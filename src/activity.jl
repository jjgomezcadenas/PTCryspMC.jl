# Activity model: the TRUE annihilation time of each event, sampled from the source's decay
# curve. (Toy / Path B: a single isotope — ¹⁵O by default. The real per-isotope model, weighted
# by the scenario budget, comes with Step 1.) This is the time the randoms pass needs; the
# detector's time-resolution blur (σ_t) is a SEPARATE step in detector.jl, not here.
#
# Decays follow the activity A(t) ∝ e^{-λ(t-t0)} over the acquisition window [t0, t1], so an
# event's annihilation time is drawn from a TRUNCATED EXPONENTIAL on [t0, t1] (rate λ =
# ln2/T½) — front-loaded, the origin of the randoms' front-loading.

"¹⁵O half-life [s] (2.037 min) — the default toy isotope."
const O15_HALFLIFE_S = 2.037 * 60.0

"""
    Isotope(name, half_life_s, beta_plus)

A positron emitter: its half-life and its positron branching fraction `β⁺` (the fraction of nuclear
decays that emit a positron, hence an annihilation). A source's activity is a NUCLEAR decay rate;
`β⁺` converts it to the annihilation rate the engine simulates.
"""
struct Isotope
    name::String
    half_life_s::Float64
    beta_plus::Float64
end

# Common PET emitters: (half-life [s], β⁺ branching). Half-lives from the usual tables.
const ISOTOPES = Dict{String,Isotope}(
    "F18"  => Isotope("F18",  109.771 * 60, 0.9686),   # clinical default
    "O15"  => Isotope("O15",  O15_HALFLIFE_S, 0.9990),
    "C11"  => Isotope("C11",  20.364 * 60,  0.9975),
    "N13"  => Isotope("N13",   9.965 * 60,  0.9981),
    "Ga68" => Isotope("Ga68", 67.71 * 60,   0.8914),
)

"Look up a PET isotope by name (e.g. `\"F18\"`); errors with the known list if absent."
function isotope(name::AbstractString)::Isotope
    haskey(ISOTOPES, name) ||
        error("unknown isotope '$name' (known: $(join(sort!(collect(keys(ISOTOPES))), ", ")))")
    ISOTOPES[name]
end

"""
    ActivityModel(; t0=0.0, t1=600.0, half_life_s=O15_HALFLIFE_S, seed=1234)

A single-isotope toy activity over the acquisition window `[t0, t1]` seconds. `event_time`
draws each event's annihilation time ∝ e^{-λ(t-t0)} on that window. `seed` keys the
deterministic per-event draws (see `event_time`). Defaults: ¹⁵O, t0 = 0 s, t1 = 600 s (10 min).
"""
struct ActivityModel
    t0::Float64        # acquisition start [s]
    t1::Float64        # acquisition end [s]
    λ::Float64         # decay constant [1/s] = ln2 / half_life
    seed::UInt64       # base seed for the per-event time draws
    norm::Float64      # 1 - e^{-λ(t1-t0)}, precomputed (the truncated-exponential normalisation)
end

function ActivityModel(; t0::Real=0.0, t1::Real=600.0,
                       half_life_s::Real=O15_HALFLIFE_S, seed::Integer=1234)
    t1 > t0 || error("acquisition window needs t1 > t0 (got t0=$t0, t1=$t1)")
    half_life_s > 0 || error("half_life_s must be > 0")
    λ = log(2.0) / half_life_s
    ActivityModel(Float64(t0), Float64(t1), λ, UInt64(seed), -expm1(-λ * (t1 - t0)))
end

"""
Build the activity model from a run config. In clinic mode (`[source].mode = \"clinic\"`) the
half-life comes from the named isotope and the acquisition window from `[source]`; otherwise it is
the toy `[timing]` activity (the count-driven phantom / API default).
"""
function ActivityModel(cfg::AbstractDict)
    if String(cfg_get(cfg, "source", "mode", "")) == "clinic"
        iso = isotope(String(cfg_get(cfg, "source", "isotope", "F18")))
        return ActivityModel(;
            t0          = Float64(cfg_get(cfg, "source", "t0_s", 0.0)),
            t1          = Float64(cfg_get(cfg, "source", "t1_s", 600.0)),
            half_life_s = iso.half_life_s,
            seed        = Int(cfg_get(cfg, "source", "time_seed", 1234)))
    end
    ActivityModel(;
        t0          = Float64(cfg_get(cfg, "timing", "t0_s", 0.0)),
        t1          = Float64(cfg_get(cfg, "timing", "t1_s", 600.0)),
        half_life_s = Float64(cfg_get(cfg, "timing", "half_life_s", O15_HALFLIFE_S)),
        seed        = Int(cfg_get(cfg, "timing", "time_seed", 1234)))
end

"""
    event_time(m, ev) -> Float64

The annihilation time [s] of event `ev`, sampled ∝ e^{-λ(t-t0)} on `[t0, t1]` (truncated
exponential, by CDF inversion). Deterministic per `(m.seed, ev)`, keyed by `event_number` — so
it is independent of how the singles were chunked across threads, identical in every pass that
reads it, and allocation-free (a stateless splitmix64 finalizer of `(seed, ev)` → a uniform,
no RNG object).
"""
@inline function event_time(m::ActivityModel, ev::Integer)::Float64
    # splitmix64 finalizer of the event id mixed with the seed → a well-distributed UInt64.
    z = UInt64(ev) + m.seed * 0x9e3779b97f4a7c15
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    z = z ⊻ (z >> 31)
    u = (z >> 11) * (1.0 / 9007199254740992.0)        # top 53 bits → uniform in [0, 1)
    m.t0 - log1p(-u * m.norm) / m.λ                   # inverse truncated-exponential CDF
end

# ── Mizuno three-component biological washout (brain) ─────────────────────────────────────────
# See CryspBrainSim/latex/washout_brain.tex. W(t) = Σ_k M_k e^{-λ_k t}, one W for all isotopes
# (tissue-level). We do NOT apply washout here (it is left to downstream); these give the per-isotope
# window-integrated survival factor g_i, stamped into each shard so downstream can apply the exact
# Bernoulli-g_i keep or recompute for a different model. Rabbit-brain params from Mizuno et al. 2003.
const MIZUNO_FRACTIONS = (0.35, 0.30, 0.35)            # M_k: fast, medium, slow
const MIZUNO_THALF_S   = (2.0, 140.0, 10191.0)         # T_k^bio [s]

# Φ(a) = (e^{a·t_irr} − 1)(e^{-a·t1} − e^{-a·t2}) / a²  — the uniform-production window integral
# (washout_brain.tex Eq. 7); a = a decay rate [1/s], the acquisition window is [t1, t2].
_washout_phi(a::Float64, t_irr::Float64, t1::Float64, t2::Float64)::Float64 =
    (exp(a * t_irr) - 1.0) * (exp(-a * t1) - exp(-a * t2)) / (a * a)

"""
    washout_gi(lambda_phys, t_irr, t1, t2; fractions=MIZUNO_FRACTIONS, thalf=MIZUNO_THALF_S) -> Float64

The per-isotope washout survival factor g_i (washout_brain.tex Eq. 8): the fraction of physically
recorded decays of a species with physical rate `lambda_phys` [1/s] that also survive clearance,
window-integrated over the acquisition [t1, t2] with uniform production over [0, t_irr]. A pure
scalar; downstream applies it as a Bernoulli keep (exact for uniform clearance).
"""
function washout_gi(lambda_phys::Real, t_irr::Real, t1::Real, t2::Real;
                    fractions::NTuple{3,Float64}=MIZUNO_FRACTIONS,
                    thalf::NTuple{3,Float64}=MIZUNO_THALF_S)::Float64
    lp = Float64(lambda_phys); ti = Float64(t_irr); a = Float64(t1); b = Float64(t2)
    denom = _washout_phi(lp, ti, a, b)
    denom == 0.0 && return 1.0
    num = 0.0
    for k in 1:3
        num += fractions[k] * _washout_phi(lp + log(2.0) / thalf[k], ti, a, b)
    end
    num / denom
end
