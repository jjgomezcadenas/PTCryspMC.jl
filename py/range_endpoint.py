"""
range_endpoint.py

Downstream analysis for an MC-vs-MC range-verification geometry study
(CRYSP open dual-head / CRYSP closed 1 m / small head PET).

The three functions below cover the parts that are independent of your
reconstruction engine. The reconstruction (identical MLEM for every arm)
sits BETWEEN thin_lm() and fit_endpoint(): you thin, reconstruct each
thinned list, form the 1-D profile, then fit its endpoint.

Pipeline per configuration:
    master LM list  --thin_lm-->  N realisations at target counts
        each realisation --> [YOUR MLEM] --> image --> 1-D depth profile
            profile --fit_endpoint--> endpoint z50 (+ z80, z20)
                {endpoints} --sigma_R--> (mean offset, sigma_R) at that dose

British spelling retained in comments for consistency with the wider project.
"""

from __future__ import annotations
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erfc


# ----------------------------------------------------------------------
# 1. Realisations by Poisson-thinning a single high-statistics master run
# ----------------------------------------------------------------------
def thin_lm(events, target_counts, n_realisations, rng=None):
    """
    True Poisson thinning: for each realisation, keep every master event
    independently with probability p = target_counts / M. The realisation
    size therefore FLUCTUATES around `target_counts` (Binomial(M, p) ~
    Poisson(target) for p << 1) — deliberately: a real acquisition at a
    given dose does not fix its count in advance, and that count noise is
    part of the sigma_R being measured. Do NOT replace this with a
    fixed-count draw; that suppresses the count fluctuation and biases
    sigma_R optimistic.

    Realisations share the one master, so they are only approximately
    independent; the approximation is good while target_counts << M.
    Rule of thumb: keep M >= 10x the largest target_counts on the dose
    axis, else sigma_R is underestimated at the high-dose end.

    Parameters
    ----------
    events : array-like, shape (M, ...)
        Master list-mode coincidences. Row = one coincidence (whatever
        columns your reconstructor needs: LOR endpoints, TOF bin, etc.).
        Only the first axis (event index) is touched here.
    target_counts : int
        EXPECTED detected-coincidence count per realisation, i.e. the
        physical count level corresponding to a chosen delivered dose.
        Map dose -> counts via: activation_yield_per_Gy * dose *
        decay/washout window fraction * geometry sensitivity. Anchor it;
        don't pick arbitrarily.
    n_realisations : int
        How many realisations to draw at this count level.
        Governs how well you estimate sigma_R (std of an N-sample std
        has relative error ~ 1/sqrt(2(N-1)); N >= 50 is a sane floor,
        a few hundred if you want tight error bars on the curve).
    rng : np.random.Generator, optional

    Returns
    -------
    list[np.ndarray]
        `n_realisations` arrays, each a row-subset of `events`, of
        VARYING length (mean target_counts, std ~ sqrt(target_counts)).
    """
    events = np.asarray(events)
    M = events.shape[0]
    if target_counts > M:
        raise ValueError(
            f"expected counts target_counts={target_counts} exceed master "
            f"statistics M={M}. Run the master with more protons, or lower "
            "the dose point."
        )
    rng = np.random.default_rng() if rng is None else rng
    p = target_counts / M
    return [events[rng.random(M) < p] for _ in range(n_realisations)]


# ----------------------------------------------------------------------
# 2. Distal-endpoint extraction by fitting the falloff (not raw R50)
# ----------------------------------------------------------------------
def _erfc_model(z, base, amp, z0, w):
    """
    Monotonic distal falloff: plateau `amp` above `base`, dropping through
    its half-height at z0 with 1/e-ish width `w`. 0.5*erfc(0)=0.5, so z0 IS
    the R50 point of the fitted (noise-averaged) edge by construction.
    """
    return base + amp * 0.5 * erfc((z - z0) / (np.sqrt(2.0) * w))


def fit_endpoint(z, profile, levels=(0.5, 0.8, 0.2), p0=None):
    """
    Fit the distal falloff of a 1-D depth-activity profile and read the
    R-level points from the fit rather than interpolating a noisy curve.
    This is the single most important choice for sigma_R in the sparse
    regime: it noise-averages the edge and yields a per-realisation
    uncertainty for free.

    Assumes `profile` runs proximal->distal in increasing z and the distal
    edge is the LAST falling edge. If your beam axis is oriented the other
    way, flip before calling.

    Parameters
    ----------
    z : 1-D array
        Depth coordinate along the beam axis (mm), monotonic increasing.
    profile : 1-D array
        Activity integrated over the fixed transverse ROI at each z.
    levels : tuple of float
        Falloff fractions of the fitted plateau to report. 0.5 primary;
        0.8 and 0.2 as secondary tail diagnostics.
    p0 : tuple, optional
        Initial (base, amp, z0, w). Auto-estimated if None.

    Returns
    -------
    dict
        'R'      : {level: z_position_mm}
        'z0'     : fitted half-height position (== R[0.5])
        'w'      : fitted edge width (mm)
        'z0_err' : 1-sigma fit uncertainty on z0 (mm)
        'popt', 'pcov' : raw fit outputs
    Returns endpoints as np.nan if the fit fails (caller should count and
    exclude failures; a high failure rate at a dose point is itself a
    result about that geometry's usability there).
    """
    z = np.asarray(z, float)
    profile = np.asarray(profile, float)

    if p0 is None:
        base0 = np.median(profile[-max(3, len(profile) // 20):])  # distal tail
        plateau0 = np.median(np.sort(profile)[-max(3, len(profile) // 10):])
        amp0 = max(plateau0 - base0, 1e-9)
        # crude half-height crossing for z0 seed
        half = base0 + 0.5 * amp0
        above = np.where(profile >= half)[0]
        z0_0 = z[above[-1]] if len(above) else z[len(z) // 2]
        w0 = 0.05 * (z[-1] - z[0]) + 1e-6
        p0 = (base0, amp0, z0_0, w0)

    try:
        popt, pcov = curve_fit(_erfc_model, z, profile, p0=p0, maxfev=20000)
        base, amp, z0, w = popt
        perr = np.sqrt(np.clip(np.diag(pcov), 0, np.inf))
        # Invert the fitted model for each requested level.
        # profile_level = base + level*amp  ->  erfc(arg)=2*level -> arg
        from scipy.special import erfcinv
        R = {}
        for lv in levels:
            arg = erfcinv(2.0 * lv)
            R[lv] = z0 + np.sqrt(2.0) * w * arg
        return dict(R=R, z0=z0, w=abs(w), z0_err=perr[2],
                    popt=popt, pcov=pcov)
    except Exception:
        return dict(R={lv: np.nan for lv in levels}, z0=np.nan, w=np.nan,
                    z0_err=np.nan, popt=None, pcov=None)


# ----------------------------------------------------------------------
# 3. Aggregate realisations -> (systematic offset, statistical sigma_R)
# ----------------------------------------------------------------------
def sigma_R(endpoints, dose_bragg_peak=None):
    """
    Collapse the per-realisation endpoints at ONE dose level into the two
    numbers that matter: the mean (systematic activity-edge position) and
    the spread (statistical precision, your discriminator).

    Parameters
    ----------
    endpoints : array-like of float
        R50 (or chosen level) from each realisation. NaNs (failed fits)
        are dropped and counted.
    dose_bragg_peak : float, optional
        Bragg-peak depth (or dose-R80) from the master run, same axis/units.
        If given, the mean endpoint is also returned as an offset from it
        (the activity-edge-to-range proxy distance; largely common-mode
        across geometries, so report it once).

    Returns
    -------
    dict with:
        n_ok, n_fail
        mean, sigma       : mean endpoint and its std (== sigma_R)
        sem               : standard error on the mean
        offset            : mean - dose_bragg_peak  (None if not supplied)
    """
    e = np.asarray(endpoints, float)
    ok = e[np.isfinite(e)]
    n_ok, n_fail = ok.size, e.size - ok.size
    if n_ok < 2:
        return dict(n_ok=n_ok, n_fail=n_fail, mean=np.nan, sigma=np.nan,
                    sem=np.nan, offset=None)
    mean, sigma = float(np.mean(ok)), float(np.std(ok, ddof=1))
    return dict(
        n_ok=n_ok, n_fail=n_fail, mean=mean, sigma=sigma,
        sem=sigma / np.sqrt(n_ok),
        offset=(None if dose_bragg_peak is None else mean - dose_bragg_peak),
    )


# ----------------------------------------------------------------------
# Sketch of the full sweep (pseudocode — plug in your reconstructor)
# ----------------------------------------------------------------------
if __name__ == "__main__":
    """
    for config in (open_dual_head, closed_1m, small_head):
        master = load_master_lm(config)            # one high-stat run
        for dose in dose_grid:
            counts = dose_to_counts(dose, config)  # physically anchored
            reals  = thin_lm(master, counts, n_realisations=200)
            eps = []
            for r in reals:
                img  = MLEM(r, config, iters=FIXED, voxel=FIXED_GRID)  # identical!
                z, prof = depth_profile(img, roi=FIXED_ROI)            # identical!
                eps.append(fit_endpoint(z, prof)['R'][0.5])
            res = sigma_R(eps, dose_bragg_peak=config.bragg_peak_depth)
            record(config, dose, counts, res)      # -> sigma_R-vs-dose curve
    # Then: overlay the three sigma_R-vs-dose curves. Closed-1m should be
    # the sqrt(N)-limited best case; dual-head the softest, and worse still
    # when the missing angle projects onto the beam axis (sweep orientation).
    """
    # Minimal self-tests.
    # (a) thin_lm is genuinely Poisson: realisation sizes fluctuate with
    #     mean = target and std ~ sqrt(target*(1-p)) (Binomial; ~ Poisson for p << 1).
    rng = np.random.default_rng(0)
    M, target = 200_000, 5_000
    sizes = np.array([len(r) for r in thin_lm(np.arange(M), target, 500, rng=rng)])
    print(f"thin_lm sizes: mean = {sizes.mean():.1f} (target {target}) | "
          f"std = {sizes.std(ddof=1):.1f} "
          f"(expected {np.sqrt(target * (1 - target / M)):.1f})")

    # (b) endpoint fit on a synthetic erfc edge with Poisson noise:
    z = np.linspace(0, 200, 400)
    truth = _erfc_model(z, base=2.0, amp=100.0, z0=150.0, w=3.0)
    noisy = rng.poisson(truth).astype(float)
    out = fit_endpoint(z, noisy)
    print(f"true R50 = 150.0 mm | fitted z0 = {out['z0']:.2f} "
          f"+/- {out['z0_err']:.2f} mm")
    print(f"R80 = {out['R'][0.8]:.2f} mm | R20 = {out['R'][0.2]:.2f} mm | "
          f"edge width w = {out['w']:.2f} mm")
