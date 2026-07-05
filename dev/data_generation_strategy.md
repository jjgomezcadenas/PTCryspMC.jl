# Data generation strategy — from proton scenario to σ_R

**Audience:** the downstream reconstruction / range-analysis instance. This explains how the
list-mode data is generated, what it means statistically, and how to consume it to produce the
range-precision figure of merit σ_R. The *directory contract* for the produced files is in
`dev/PRODUCTS.md`; the *endpoint-extraction method* is in `docs/range_verification_recipe.md`.
This document is the bridge between them.

You do **not** need the simulator (`PTCryspMC.jl`) to do your work — you consume the products tree
(`PtCryspProds/`) and the frozen scenarios. But understanding how the products were made is
essential to not misuse them (especially the shard/thinning distinction below).

---

## 1. The pipeline and the four parts

```
ptcrysp-scenarios   →   PTCryspMC.jl   →   PtCryspProds/   →   THIS repo (recon + σ_R)
 (proton+phantom       (photon MC,        (list-mode LORs,     (reconstruct, thin,
  → annihilation        detector           per scanner/         fit endpoint, σ_R)
  points, frozen)       response)          crystal, sharded)
```

- **`ptcrysp-scenarios`** — a frozen Geant4 proton run: a SOBP field on a phantom, producing the
  *spatial* positron-annihilation source (`emitters.csv`) + a per-isotope decay budget. One
  scenario = one (proton field + phantom). Detector-independent.
- **`PTCryspMC.jl`** — a photon-only Monte Carlo: draws annihilations from the scenario, transports
  the 511 keV pairs through the phantom + a scanner, applies the detector response, and writes the
  coincidence (LOR) list `lors_det.h5`. Detector-specific.
- **`PtCryspProds/`** — the produced LOR files, organized by scenario / scanner / crystal /
  timing-dose, stored as **shards** (see §3). This is your input.
- **This repo** — reconstructs each realization, extracts the distal range endpoint, and reports
  σ_R vs dose, one curve per scanner geometry.

---

## 2. What one LOR file is

Each `lors_shardNNN.h5` (schema: `PTCryspMC.jl/docs/SCHEMA.md`) is the list of coincidences **one
scanner+crystal would record from one SOBP irradiation** — per accepted coincidence: two 3-D hit
points (mm, detector-smeared), two energies (keV), two times (ns), per-gamma phantom-scatter
counts, a timing residual, and a **truth flag** (0 = true, 1 = scatter, 2 = random). Root
attributes carry the full provenance (scenario, scanner, crystal, dose, budget, master_seed, shard
index, detector windows, geometry segmentation).

The detector response is already baked in: position resolution (σ_xyz), energy resolution and the
photopeak cut, the coincidence-time window. You receive *detected, selected* LORs — you do not
re-apply the detector model.

---

## 3. THE CENTRAL CONCEPT: shards vs realizations (two sampling layers)

This is the one thing you must get right. There are **two** independent random layers, and they
must not be conflated.

### Layer 1 — shards (upstream, already done, stored on disk)

The simulator produces, per (scenario, scanner, crystal, budget, top-dose), a set of **shards**:

- A shard = one **full independent MC run** — its own source Poisson draw (`M_j ~ Poisson(N_j)`
  annihilation points) **and** its own photon transport (emission directions, interactions).
- There are **N_shard ≈ 10** of them (`lors_shard000.h5 … lors_shard009.h5`), each ~10⁸ decays at
  the **top dose** of the study.
- **Pooled, the 10 shards are the master** (~10⁹ decays). Statistically identical to one monolithic
  10⁹ run (all iid samples from the emitter pool), but built as 10 production-scale runs to keep
  each run's memory/time bounded.
- **10× rule:** the master must be ≥ 10× the largest count you will thin to, so the top dose point
  thins with keep-probability p ≤ 0.1. Ten shards at the top dose give exactly that headroom.

A shard is seeded upstream by `(master_seed, shard_index)`. **The same shard index across scanners
and crystals sees the identical source** (the annihilation points depend only on the scenario +
seeds, never on the detector). That common-mode is what makes the geometry and detector comparisons
fair — see §6.

### Layer 2 — realizations (downstream, YOUR job, never stored)

A **realization** is one simulated "experiment" at a chosen dose, produced by **Poisson-thinning the
pooled master**. You generate Z of them (Z ≈ 100–200) *in memory*, reconstruct and fit each, and
discard it. **Realizations are never written to disk.** σ_R = the spread of the fitted endpoints
across the Z realizations.

**Do not** treat a shard as a realization for the main curve, and **do not** call your thinned
subsamples "shards." The stored files are shards; your transient draws are realizations.

---

## 4. `thin_lm` — the thinning you must implement (a few lines)

Write it fresh (do not port any exact-count version — that contradicts this design). It is a seeded
Bernoulli stream filter over the **pooled** shards:

```
thin_lm(shard_files, target_counts, realization_index):
    M_total = Σ nrows(f) for f in shard_files          # the pooled master size
    p       = target_counts / M_total                   # keep-probability (≤ 0.1 at the top dose)
    rng     = seeded by realization_index               # OWN namespace (not the upstream seeds)
    for each event e streamed across ALL shard_files:
        if rng.uniform() < p:  yield e                  # keep independently, per event
```

Three non-negotiable properties (each fixes a specific failure mode):

1. **Pool across all shards.** `p = target/M_total` where `M_total` sums over *all* shard files, and
   the mask spans the union. Thinning each shard file separately is the "p=1 per file" degenerate
   trap — it would return whole shards, not a realization.
2. **Bernoulli, not exact-count.** Keep each event independently with probability p; the realization
   count then **fluctuates** as Binomial(M_total, p) ≈ Poisson(target). That count fluctuation is
   part of σ_R — a fixed-count draw suppresses it and biases σ_R optimistic. Do not force an exact
   count.
3. **Own seed namespace.** Seed by the *downstream* realization index only. It has nothing to do
   with the upstream `(master_seed, shard_index)`; it operates purely on the produced LORs.

**Thin down, never up.** The shards are at the master's top dose; you produce realizations at any
dose ≤ that by choosing `target_counts = dose_to_counts(dose)`. p stays ≤ 0.1 by the 10× rule. There
is no "master dose" to exceed.

Anchor `dose_to_counts` physically (activation yield/Gy × dose × decay/washout window fraction ×
geometry sensitivity); see `docs/range_verification_recipe.md` step 2.

---

## 5. The σ_R pipeline (per configuration)

For one (scenario, scanner, crystal, budget):

```
pool = glob(<...>/lors_shard*.h5)                    # the master for this configuration
for dose in dose_grid:
    counts = dose_to_counts(dose)
    eps = []
    for z in 1..Z:                                    # Z realizations, in memory
        lm   = thin_lm(pool, counts, realization_index=z)
        img  = MLEM(lm, scanner_geometry, mu_map, iters=FIXED, voxel=FIXED)   # identical across arms
        z_prof = depth_profile(img, roi=FIXED)        # 1-D activity(z) over a fixed transverse ROI
        eps.append(fit_endpoint(z_prof)['R'][0.5])    # erfc-fit R50 (range_endpoint.py)
    record(dose, sigma_R = std(eps), mean = mean(eps))
```

- **MLEM** and **depth_profile** are yours to write. **`fit_endpoint`/`sigma_R`** already exist in
  `range_endpoint.py` (migrating into this repo). Reconstruction inputs beyond the LORs: the
  **scanner geometry** (system model) and the **phantom μ-map** (attenuation correction) — both
  shipped in the products tree (`PtCryspProds/PRODUCTS.md`).
- **Fix everything except geometry** across arms: same voxel grid, same ROI, same MLEM iteration
  count, same endpoint code. These live in one `config/` and are consumed unchanged — a per-arm
  inconsistency here silently corrupts the comparison.
- **Trues-only on the first pass** (truth flag == 0). Add scatter/randoms as a separate axis later;
  open and closed geometries handle scatter differently and it must not confound the headline curve.

---

## 6. The two comparisons (why the common-mode source matters)

- **Geometry comparison (the headline):** same scenario / crystal / budget, sweep the **scanner**
  directories. One σ_R-vs-dose curve per scanner. Because matched shard indices carry the identical
  source, the arms differ only by the scanner — the cross-section and source systematics are
  common-mode and cancel (this is the MC-vs-MC argument that lets us use the simplest endpoint
  metric; `docs/range_verification_recipe.md`).
- **Detector comparison:** same scanner / budget, sweep the **crystal** directories.

To keep the common-mode exact, thin the *same* realization index across arms from *matched* shards.

---

## 7. The independent-shards cross-check

The 10 shards are themselves 10 truly independent experiments at the top dose. Fit each shard's
endpoint directly and take the std of the 10 — a **bias-free** σ_R at the top dose (no shared-master
deflation), coarse (n=10 → ~24% uncertainty on the std). The thinned σ_R at the top dose must agree
with it within that precision. Disagreement means the thinning (or the master size / the 10× rule)
is wrong. Cheap (ten fits) and worth running as a standing gate.

---

## 8. What this repo must build vs reuse

- **Build:** `thin_lm` (§4), `MLEM` (with attenuation correction from the μ-map + the scanner system
  model), `depth_profile` (3-D image → 1-D activity(z) over a fixed ROI), the sweep driver (§5).
- **Reuse:** `range_endpoint.py` — `fit_endpoint` (erfc R50 from the falloff) and `sigma_R`
  (aggregate). Note its current `thin_lm` is a stub; replace it per §4.

---

## 9. Concrete numbers (the first available case)

- Scenario `uniform_headep_sobp_1e8`: SOBP on a uniform G4_BRAIN_ICRP head, 2.10M annihilation
  points in the pool; 5 isotopes (O15, C11, N13, C10, O14); budgets fast / inroom / offline.
- One shard at 1 Gy (fast) = **80.18M decays** → **17.4M photopeak-selected LORs** (BGO, closed
  ring), 81% true / 19% scatter / 0.18% random. Full detector response applied (σ_xyz 1.7 mm, 10%
  energy FWHM, 450 keV photopeak cut, τ = 3 ns).
- A full master = **10 such shards** (~50 min/detector, or run in parallel). Currently **one shard
  exists** (shard 0, BGO). The rest (shards 1–9, the CsI arm, other scanners) are additional
  production runs into the same tree — no new machinery.

---

## 10. Gates (arms not yet producible)

Two scanner arms exist as directory slots but cannot be generated until the simulator gains the
matching capability:
- **Open dual-head** — needs partial-ring (`phi_gaps`) / planar-panel geometry in the engine.
- **Mixed crystal** (e.g. BGO core + CsI wings) — needs per-block crystal materials in the engine.
Closed-ring, head-sized, and children's scanners are `CylShell` variants and are producible now
(config-only). Do not expect the gated arms in `PtCryspProds/` yet.
