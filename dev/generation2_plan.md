# Generation-2 plan — centred phantom, irradiation-referenced timing, isotope truth

Products generation tag: **`generation = "v2"`** (stamped in every shard; see §5).

**Status:** plan (2026-07-12). Supersedes the off-centre survey masters (kept as legacy — see
`dev/status.md` "Source centring"). This document is the **downstream contract** for what the new
list-mode products carry and mean. Read it with `dev/data_generation_strategy.md` (the shard/
realization/σ_R model, unchanged) and `dev/PRODUCTS.md` (the directory layout).

Three things are fixed in this generation, all in this repo, so downstream receives clean products:

1. **Phantom position** — the patient is placed so the activity distal edge (range endpoint, R50)
   sits at the ring centre z = 0, not off-centre.
2. **Timing** — the acquisition is a fixed-length window sliding with the patient-arrival delay,
   on an **irradiation-end clock**, delivered as a family of fixed scenarios.
3. **Downstream information** — each LOR carries its **emitting isotope** (plus the true origin and
   the decay time already present), so downstream can do exact per-isotope work and apply washout
   itself. Strictly documented here and in `docs/SCHEMA.md`.

---

## 1. Fix 1 — phantom centring (tumour)

The frozen scenario puts the tumour off the ring centre (native frame: SOBP-target centre z ≈
−25 mm). Clinically the patient is positioned with the **tumour** at isocentre — a fixed anatomical
reference the clinician sets from NMR/CT *before* choosing any acquisition window — so we rigidly
shift the **emitters and the phantom together** in z (the scanner stays at the origin) until the
dose-target centre sits at z = 0.

- Knob: `[source].center_on = "tumour"` (`src/scenario.jl`), the v2 default. The target centre is the
  distal dose R80 (the clean range endpoint in z) minus half the target thickness `(dist−prox)/2`
  from run_meta — dose-based, so it is **independent of the acquisition window / isotope mix /
  washout**: one fixed offset (≈ +25.6 mm here) for every scanner, crystal, and scenario. Stamped
  `source_z_offset_mm` + `center_on`.
- With the tumour at isocentre the activity R50 lands ~+10 mm and the dose R80 ~+20 mm distal of
  centre — negligible for acceptance in these AFOVs.
- `center_on = "distal_edge"` (centre the activity R50) remains in the code but is NOT the default: it
  drifts with the window's isotope mix, so it is not a stable positioning reference.
- All v2 products are tumour-centred. Legacy off-centre masters keep the old `<budget>_<dose>` leaf
  token; v2 products use the timing token (§2.5), so the two never collide.

## 2. Fix 2 — irradiation-referenced timing

### 2.1 The clock

`t = 0` is the **end of irradiation** (beam off). Every decay time is measured from there. The LOR
column `t_decay_s` therefore has **zero = irradiation end** (was: acquisition start). One physical
consequence that removes an entire class of downstream error: on this clock the acquisition window
is a pure band, washout is a pure per-event function of absolute time, and physical decay is in the
event density — the three factor and do not couple.

### 2.2 The scenarios (fixed, clean deliverables)

The patient arrives at the scanner a delay `t_del` after the beam ends and is imaged for a fixed
acquisition `t_ac`. Each scenario is the window **[t1, t2] = [t_del, t_del + t_ac]** — a
**full-length `t_ac` acquisition** whose start slides with `t_del` (not a shrinking window).

| scenario | t_del [s] | t_ac [s] | window [s] |
|---|---|---|---|
| del120 | 120 | 300 | [120, 420] |
| del180 | 180 | 300 | [180, 480] |
| del300 | 300 | 300 | [300, 600] |

(t_ac = 300 s ≈ realistic in-room PET; 120–180 s realistic patient-move, 300 s late. t_ac and the
set are config-driven and extensible.)

### 2.3 Transport once, band-cut per scenario

Transport is the only expensive step, so it runs **once** over the **union window** [120, 600] s
(min t_del to max t_del + t_ac);
the per-scenario shards are cheap band-cuts of that single transport (the decay time is assigned in
the coincidence builders, so the transport itself is window-agnostic):

1. Sample `M_j ~ Poisson(N_j(union) · f_inside_j)` and transport once → union singles + `lors_truth`
   + `randoms`, each event carrying its decay time on the irradiation-end clock and its isotope.
2. For each scenario, **band-cut** `t_decay ∈ [t_del, t_del + t_ac]`, smear + energy-select, and
   write one `lors_det.h5` per scenario. No re-transport.

Band-cutting the union transport is statistically identical to a per-scenario transport (correct
counts and Poisson noise per window); it just avoids ~3× the transport cost.

### 2.4 The counts: N_j⁰ and the derivation

The upstream budgets (`fast`/`inroom`/`offline`) are window-integrals of one decaying
irradiation-end population `N_j⁰`. We back `N_j⁰` out of the `fast` budget analytically (production
ends before any window, so pure decay):

```
N_j⁰      = N_j(fast) / (e^{-λ_j · t_del_fast} − e^{-λ_j · (t_del_fast + t_meas_fast)})   # [60,1260]
N_j(win)  = N_j⁰ · (e^{-λ_j · t1} − e^{-λ_j · t2})                                          # any [t1,t2]
```

Cross-check (done): `N_j⁰` from `fast` and from `inroom` agree to < 0.2 % per isotope (O15 7.80e7,
C11 4.01e7, C10 1.57e6), confirming the budgets are one population — so **no upstream round-trip is
needed** to generate any (t_del, t_ac). `λ_j = ln2 / T½` from `isotopes.csv`.

### 2.5 Products layout (per scenario)

The timing scenario replaces the budget token in the leaf:

```
<scenario>/<scanner>/<crystal>/del<t_del>s_ac<t_ac>s_<dose>/lors_shard<NNN>.h5
e.g. uniform_headep_sobp_1e8/crysp_ring_1m_csi_2x0/csi/del60s_ac1200s_1Gy/lors_shard000.h5
```

Three leaves per (scanner, crystal), one per scenario, each a 10-shard master. The shard/realization/
`thin_lm`/σ_R model of `data_generation_strategy.md` is **unchanged**: it applies per scenario leaf
(pool the 10 shards of one leaf, thin down, N× rule). `truth/`, `phantom/`, `scanner_geometry.json`
are shared as before.

## 3. Fix 3 — isotope truth (what downstream carries on with)

Each LOR now carries the emitting **isotope id** (one per LOR; both gammas share the annihilation).
Together with fields already present, every LOR gives:

| field | meaning |
|---|---|
| `isotope` | emitting isotope id (0=O15, 1=C11, 2=N13, 3=C10, 4=O14) — **new** |
| `x0/y0/z0_mm` | true annihilation point (origin depth z0) — already present |
| `t_decay_s` | decay time, **zero = irradiation end** — semantics changed (§2.1) |
| `truth` | 0 true / 1 scatter / 2 random — already present |

This is strictly more than the label-free reconstruction can recover. It enables:
- **Exact** per-species σ_R (pure sub-samples, not the enriched lower bound of an isotope-blind
  posterior);
- **validation** of the label-free posterior `P(i | z0, t_d)` against truth;
- **blind-detector emulation** whenever wanted — just ignore the isotope column.

The isotope column propagates singles → `lors_truth` → `lors_det` and into `randoms`; documented in
`docs/SCHEMA.md` (regenerated, a test guards drift).

## 4. Washout — left to downstream (per `CryspBrainSim/latex/washout_brain.tex`)

We do **not** thin for washout. The note proves a per-species Bernoulli keep with `g_i` is exact
(mean profile *and* Poisson σ_R), position-independent, and equally valid applied downstream; and
because the isotope label is now carried, downstream applies the exact `g_i` keep directly (its
recommended path), able to vary the large-uncertainty Mizuno parameters for free. We ship the
un-washed physical-window counts plus everything needed to apply washout:

- **Model:** Mizuno three-component brain washout, one `W(t) = Σ_k M_k e^{−λ_k t}` for all isotopes
  (tissue-level): fast M=0.35 T=2.0 s, medium M=0.30 T=140 s, slow M=0.35 T=10191 s.
- **Stamped** in each shard as provenance/convenience: `t_irr_s`, `t_del_s`, `t_ac_s`, the Mizuno
  `(M_k, T_k)`, and the per-isotope closed-form `g_i` for the shard's window (note Eq. 8:
  `g_i = Σ_k M_k Φ(λ_i+λ_k)/Φ(λ_i)`, `Φ(a)=(e^{a·t_irr}−1)(e^{−a·t1}−e^{−a·t2})/a²`).
- Downstream can also apply the age-resolved `W(τ)` per event (τ = `t_decay_s`, the age since
  irradiation end) — statistically identical for uniform W, and the hook for a future spatial
  `W(τ; x)`. The un-washed per-species profiles `P_i(z)` it needs live in the physical `truth/`
  bundle already shipped.

## 5. Self-describing metadata (root attributes)

Every `lors_det.h5` carries the full context in its HDF5 root attributes, so a shard is interpretable
standalone — no external doc needed to read it. The drivers set them; `copy_provenance!` carries the
propagated ones down singles → `lors_truth` → `lors_det`; and the attribute set is **enumerated in
`docs/SCHEMA.md` and covered by the same drift test as the columns** (a missing/renamed attr fails a
test). Types: `str`, numeric scalar, or array. **`generation="v2"` is mandatory** — consumers must
check it and refuse to mix generations, since the `t_decay_s` zero moved (§2.1).

**A. Generation & schema** (the semantics guard):
`generation="v2"`, `schema_version` (int, drift-guarded), `t_decay_zero="irradiation_end"`,
`has_randoms` (bool), `washout_applied=false` (bool).

**B. Source & placement:**
`source_mode="api"`, `scenario`, `budget_source="fast"` (the budget N_j⁰ was derived from),
`dose_Gy`, `master_seed`, `realization` (shard index), `center_on="distal_edge"`,
`source_z_offset_mm`, `isotope_names` (str[]), `isotope_half_lives_s` (f64[]),
`n_expected` (f64[]: per-isotope mean decays in this scenario window).

**C. Timing** (irradiation-end clock):
`t_irr_s`, `t_del_s`, `t_ac_s`, `t1_s`, `t2_s` (this scenario window = [t_del, t_del+t_ac]),
`t_union_lo_s`, `t_union_hi_s` (the transport window), `time_seed`.

**D. Detector & selection:**
`crystal`, `reco_emin_keV`, `emin_keV`, `tau_ns`, `eres_fwhm_511`, `sigma_xyz_mm`, `sigma_t_ns`,
`prompt_gamma_modeled`.

**E. Geometry (full self-containment):**
`geometry_file` (name), `scanner_name`, `ring_radius_mm`, `afov_mm`, `crystal_depth_mm`, `n_phi`,
`n_z`, and **`geometry_json`** — the entire geometry JSON (world + phantom + scanner) embedded
verbatim as a string, so one `lors_det.h5` alone fully reconstructs its scanner. The
`scanner_geometry.json` sidecar remains the human-readable copy; the embedded string is authoritative
for the shard.

**F. Washout inputs** (for downstream to apply `g_i`; **not** applied here):
`washout_model="mizuno_brain_3comp"`, `washout_fractions` (f64[] = M_k), `washout_Thalf_s`
(f64[] = T_k), `washout_g` (f64[]: per-isotope g_i for this window). `washout_g` is **computed from**
the stamped Mizuno params + timing (note Eq. 8), never an independent number, so the file cannot
self-contradict; downstream may recompute or apply the age-resolved W(τ).

**G. Counts:** `nevents` (transported annihilations), `nrows` (LORs in this file).

Large data (truth profiles, μ-map) stays in the shared `truth/`, `phantom/` bundles — attributes are
metadata (scalars, small arrays, short strings, the geometry JSON), never bulk datasets.

## 6. What changed vs the off-centre generation

| | off-centre (legacy) | generation-2 |
|---|---|---|
| phantom | native off-centre | distal edge at z=0 (`center_on`) |
| `t_decay_s` zero | acquisition start | **irradiation end** |
| timing | one budget window | fixed scenario family (t_del × t_ac) |
| leaf token | `<budget>_<dose>` | `del<t_del>s_ac<t_ac>s_<dose>` |
| isotope in LOR | no | **yes** |
| washout inputs stamped | no | Mizuno + `g_i` + timing |

Legacy masters are kept (distinct leaf token); the σ_R comparison is on generation-2.

## 7. Downstream action items

- `docs/SCHEMA.md` — regenerated: new `isotope` column; `t_decay_s` zero semantics; new attrs.
- `dev/PRODUCTS.md` — leaf naming (§2.5).
- `dev/data_generation_strategy.md` — shard/realization/`thin_lm`/σ_R **unchanged**; note it now runs
  per scenario leaf, and add the washout `g_i` step (from the stamped attrs / `washout_brain.tex`).

## 8. Alternative considered (not chosen)

Ship **one** union-window shard carrying `t_decay_s` and let downstream band-cut per scenario (¼ the
storage, band-cut trivial on the irr-end clock). Not chosen: the request is fixed-scenario clean
counts here. Each shard still carries `t_decay_s`, so downstream can refine within a scenario if it
wants.

## 9. Build order

1. `src/scenario.jl` — read `t_irr_s`, `t_del_s` from budget meta; derive `N_j⁰`; set the union
   window [120,600] + `N_j(union)`; activity models over the union window (irr-end clock).
2. `src/coincidences_hdf5.jl` — add the `isotope` LOR column (+ doc, provenance attrs).
3. build_true_coincidences / build_randoms — carry isotope; assign `t_decay` on the union window.
4. `scripts/reco_lors.jl` — loop scenarios: band-cut, emit one `lors_det.h5` per scenario; stamp the
   full §5 attribute set (`generation="v2"`, `schema_version`, `t_decay_zero`, timing, geometry incl.
   embedded `geometry_json`, Mizuno + computed `g_i`).
5. `scripts/run/publish_prod.jl` — per-scenario leaf token `del<t_del>s_ac<t_ac>s_<dose>`.
6. `scripts/gen_schema.jl` → regenerate `docs/SCHEMA.md` (columns **and** the §5 attribute set),
   extend the drift test to the attributes; `Pkg.test`.
7. First shard: `crysp_ring_1m` CsI, shard 0 → the three scenario `lors_det` (del120/180/300).
