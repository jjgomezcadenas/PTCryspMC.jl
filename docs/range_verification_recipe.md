# Distal-endpoint range-verification study — analysis recipe

**Goal.** Compare three PET geometries for proton-therapy range verification on an
MC-vs-MC basis, using the distal activity endpoint as the metric:

1. **CRYSP open** — dual-head "in-beam", one head above / one below the patient.
2. **CRYSP closed** — long "in-room" system, closed ring, AFOV ~1 m.
3. **Small head PET** — standard compact scanner, reference point.

**Terminal figure of merit.** Statistical precision on the distal activity
endpoint, `sigma_R`, as a function of delivered dose — one curve per geometry.
This is the range-verification analogue of the clinical CRC-vs-BV curve.

**Why distal endpoints (not filtering / MLS).** MC-vs-MC makes cross-section
systematics common-mode across the three arms, so they cancel in the geometry
comparison. That removes the usual objection to endpoint methods (absolute-activity
MC dependence) and lets us use the simplest, most telling metric. Filter-function
and most-likely-shift methods buy robustness we don't need here and add machinery
that obscures the geometric response we're trying to isolate.

---

## The chain: from LM LORs to sigma_R-vs-dose

### 1. Reconstruct — identically for all three geometries
The limited-angle penalty of the open dual-head arm **is a reconstruction
artefact**: it does not exist in the raw LOR list, only when an image is recovered
from incompletely sampled angles. So do **not** shortcut from LM LORs straight to a
profile for the open geometry — that erases the effect being measured. Use **MLEM**:
consistent across arms, and it degrades gracefully under missing angles where FBP
would streak (the fair way to give the open geometry its best honest shot).

### 2. One master run -> many realisations by Poisson-thinning
Do **not** run N independent MC simulations per dose level. Run one high-statistics
master per configuration, then thin the list-mode stream — each event kept
independently with probability p = target/M — many times over per target dose. Each
realisation's count fluctuates around the target, exactly as a real acquisition at
that dose would (that count noise is part of the sigma_R being measured — do not
fix the count); the whole dose axis comes from a single expensive run. Keep the
master well above the largest dose point (M >= 10x its counts), else realisations
stop being independent and sigma_R is underestimated there. Anchor the mapping
physically:

`dose -> counts  =  activation_yield_per_Gy x dose x (decay+washout window fraction) x geometry_sensitivity`

**Sharded master.** The top of the dose axis sets the master size: with a top dose
point worth 10^8 decays, the master needs 10^9 (the 10x rule above). Build it as
ten independent 10^8-decay runs with different seeds — statistically identical to
one monolithic 10^9 run, and every simulation, file and pruning step stays at the
proven production scale. A realisation is then drawn from the POOLED master: stream
all ten shards and keep each event independently with p = target/M_total (p = 0.1
at the top dose point, smaller below). Seed the thinning RNG with the realisation
index: each realisation is reproducible on demand, so it can be reconstructed,
fitted and discarded — no realisation is ever stored.

### 3. Reconstruct each realisation — everything fixed
Same MLEM iteration count, same voxel grid (identical across geometries despite
different FOVs), no post-filter on the first pass. Iteration count is common-mode
because distal-edge sharpness depends on it.

### 4. Collapse to the 1-D depth-activity profile
Integrate over a **fixed transverse ROI** around the beam axis. Same ROI definition
for all arms — this is where a per-geometry inconsistency most easily sneaks in.

### 5. Extract the endpoint by fitting the falloff (not raw R50)
Reading 50%-of-distal-max off a noisy curve inherits all the noise. Fit a
complementary-error-function (or logistic) edge to the falloff region and take R50
analytically from the fit. Gives a continuous, noise-averaged endpoint plus a
per-realisation fit uncertainty, and makes R80/R20 usable as secondary tail
diagnostics. **This single choice does more for `sigma_R` than anything downstream.**

### 6. Aggregate across realisations -> the two numbers that matter
At each dose level:
- **mean endpoint** = systematic (the activity-edge position; largely common-mode
  across geometries — report once, as the activity-edge-to-Bragg-peak offset).
- **std across realisations** = `sigma_R`, the statistical precision — the
  **discriminator**.

Sweep dose -> `sigma_R`-vs-dose curve, one per configuration.

---

## Guardrails

- **Cross-check the thinning at the top dose.** The ten shards are themselves ten
  truly independent experiments at exactly the top dose point (10^8 each): fit the
  endpoint of each shard and take the std of the ten — a bias-free sigma_R with no
  shared-master deflation. The thinned estimate at that dose must agree with it
  within the ten-sample precision (~1/sqrt(2x9) ~ 24%); agreement validates the
  thinning machinery at the one dose point where p is largest. Costs ten fits.
- **Trues-only on the first pass.** Add scatter/randoms as a separate axis later:
  open and closed geometries handle scatter qualitatively differently, and that must
  not be confounded into the headline curve.
- **Sweep beam-axis orientation for the open geometry.** Endpoint precision depends
  on whether the missing angle projects onto the beam direction; a single orientation
  will either flatter or unfairly damn the dual-head arm.
- **Fix everything except geometry** — same isotope cocktail, activation scoring,
  acquisition/decay window, reconstruction, voxel size, ROI, endpoint-extraction
  code — so `sigma_R` differences are purely geometric.

## First sanity check (before the full sweep)
Reconstruct one high-count master per geometry; overlay the three 1-D profiles with
the dose curve. If the closed-1 m distal edge is not visibly sharpest and the
dual-head not visibly softest, something in the reconstruction or ROI is off — catch
it before generating hundreds of realisations.

## Expectations
- **Closed 1 m** — full angular coverage, sqrt(N)-limited best case; the reference
  the other two are scored against.
- **Open dual-head** — softest, and worse still when the missing angle projects onto
  the beam axis (limited-angle elongation).
- **Small head PET** — tests whether compact high-solid-angle-per-cost competes with
  raw AFOV for a localised (head-sized) activation volume; may beat its size.

## Governing scaling
Range precision scales as `sigma_R ~ 1/sqrt(N_coincident_counts)`. Counts-limited,
not resolution-limited — so the long-AFOV geometric sensitivity feeds straight into
precision, and (unlike the clinical CRC case) the missing TOF hurts less because the
analysis is 1-D along a known beam axis and the distal edge is strongly localised.
The same scaling sizes the master: the top dose point dictates M >= 10x its counts
(10^8-decay top point -> 10^9-decay master, i.e. the ten shards above), while every
lower dose point thins from the same master for free.

## Companion code
`range_endpoint.py` — `thin_lm` (step 2), `fit_endpoint` (step 5), `sigma_R` (step 6).
Reconstruction + `depth_profile` (step 4) slot in between; the MLEM step is yours.
