# PtCryspProds — the products directory contract

The layout of the list-mode products produced by `PTCryspMC.jl` and consumed by the downstream
reconstruction / range-analysis repo. This is the interface between the two: the simulator writes
this tree (via `publish_prod`), the analysis repo reads it. For *why* the data is shaped this way
(shards, thinning, σ_R), read `dev/data_generation_strategy.md` first.

`PtCryspProds/` is a directory **sibling to the repos** (like `ptcrysp-scenarios/`), not committed
to `PTCryspMC.jl` (the `.h5` files are large and regenerable). It is the authoritative store of the
simulated LOR data.

---

## The tree

```
PtCryspProds/
  README.md                              self-describing copy of this contract (written by publish_prod)
  SCHEMA.md                              the shard column schema (docs/SCHEMA.md; refreshed on every publish)
  scanner_prods.pdf                      the scanner-productions note: systems, CTR calibration, deliverables
                                         (latex/scanner_prods.tex; refreshed on every publish)
  <scenario>/                            ── SCENARIO = proton run + phantom (fixes the SOURCE + phantom)
    phantom/                             ── shared by ALL scanners/crystals — the μ-map inputs
      phantom_regions.csv                   the region(s): shape, semi-axes, centre, material (from the scenario)
      material_<name>.csv                   composition + μ(511) per material (build the voxel μ-map from these)
    truth/                               ── shared by ALL scanners/crystals — detector-independent truth
      depth_dose.csv                        dose(z): the Bragg/SOBP distal edge → dose-R80
      sobp_layers.csv (+ _meta)             the SOBP beam design
      run_meta.csv                          target-box depths, Np/Gy, normalization
      sampling_budget_<budget>.csv (+_meta) per-isotope N_expected for the acquisition timing
      activity_profile_<budget>.csv (+_meta) binned true activity(z), per isotope + total → activity-R50
    <scanner>/                           ── SCANNER GEOMETRY, named incl. crystal depth in X0 when it matters
                                            (crysp_ring_1m_bgo_2x0, crysp_ring_1m_csi_2x0, …; also head,
                                            children, open_dualhead, mixed…)
      scanner_geometry.json                 the ring/panel geometry (+ per-block crystal map if heterogeneous)
      <crystal>/                         ── CRYSTAL (bgo, csi) — ONLY for homogeneous scanners
        <budget>_<dose>/                 ── acquisition timing budget + the master's TOP dose
          config.toml                       the regeneration recipe (identical across shards but shard_index)
          lors_shard000.h5                  shard 0  ┐
          lors_shard001.h5                  shard 1  │  the ~10 shards = the master
          …                                          │  (pool + thin downstream → σ_R)
          lors_shard009.h5                  shard 9  ┘
```

**Heterogeneous scanner** (crystal is intrinsic, e.g. BGO core + CsI wings): no `<crystal>/` level —
`<budget>_<dose>/` hangs directly under `<scanner>/`, and the per-block crystal map lives inside its
`scanner_geometry.json`.

Example (the first available case):
```
PtCryspProds/uniform_headep_sobp_1e8/
  phantom/{phantom_regions.csv, material_g4_brain_icrp.csv}
  truth/{depth_dose.csv, sobp_layers.csv, run_meta.csv, sampling_budget_fast.csv, activity_profile_fast.csv, …}
  crysp_ring_1m/
    scanner_geometry.json
    bgo/fast_1Gy/{config.toml, lors_shard000.h5 … lors_shard009.h5}
    csi/fast_1Gy/{config.toml, lors_shard000.h5 … }
```

---

## The axes

| level | axis | varies | shared with… |
|---|---|---|---|
| `<scenario>/` | proton field + phantom | scenario | the **source**, the **μ-map** (`phantom/`) & the **truth** (`truth/`) — across everything below |
| `<scanner>/` | scanner geometry | scanner | `scanner_geometry.json` — across the crystals in it |
| `<crystal>/` | crystal material (homogeneous only) | crystal | — |
| `<budget>_<dose>/` | acquisition timing + master top dose | budget | — |
| `lors_shardNNN.h5` | **shard index** (a master component) | shard | matched by index across scanners/crystals → identical source |

**The source is common-mode across every scanner and crystal at a matched shard index** — the
annihilation points depend only on `(scenario, master_seed, shard_index)`, never on the detector.
That is what makes the geometry and detector comparisons isolate a single axis.

---

## Definitions (read these before using the tree)

- **Scanner** = a ring/panel geometry **plus its crystal assignment**. Homogeneous: one material,
  swapped via `config.crystal_material`, and the tree has a `<crystal>/` level. Heterogeneous: a
  per-block crystal map baked into `scanner_geometry.json`, and there is **no** `<crystal>/` level.
- **Shard** (`lors_shardNNN.h5`) = one full independent MC run at the top dose; ~10 shards pooled =
  the master. A shard is **NOT** a σ_R realization — see below.
- **Realization** = a downstream, in-memory Poisson-thinned subsample of the pooled master. Produced
  by `thin_lm`, reconstructed, fitted, discarded. **Never stored here.** (The upstream
  `[source].realization` config field is the shard's source seed — labelled the *shard index* in
  this tree.)
- **Budget** = the acquisition-timing scenario (fast / inroom / offline); different budgets are
  genuinely different sources (different N_j and randoms) — you cannot thin between them.
- **Dose** in the leaf name = the master's **top** dose. Lower doses are produced downstream by
  thinning *down*; they are not stored.

---

## Reading rules (for the analysis repo)

- **The master for one configuration:** `glob(<scenario>/<scanner>/<crystal>/<budget>_<dose>/lors_shard*.h5)`
  → pool all shards → `thin_lm` (Bernoulli p = target/M_total over the union; see
  `dev/data_generation_strategy.md` §4) → realizations → σ_R.
- **Geometry comparison (headline):** fix scenario/crystal/budget, sweep `<scanner>/`.
- **Detector comparison:** fix scanner/budget, sweep `<crystal>/`.
- **Reconstruction inputs:** `<scenario>/phantom/` (build the μ-map) + `<scanner>/scanner_geometry.json`
  (the system model). Never in the LOR file.
- **Truth reference (scoring):** `<scenario>/truth/` — detector-independent. `activity_profile_<budget>.csv`
  is the clean β⁺ source curve (→ activity-R50, the recon target); `depth_dose.csv` is the physical dose
  edge (→ dose-R80). Their offset is the locked reference the reconstructed edge is scored against. The
  activity profile carries the *same* source scaling as the shards (`N_expected·f_inside` at the run's
  dose; escaped positrons excluded), so it composes with the pooled LORs directly.
- **Provenance / regeneration:** every `lors_shardNNN.h5` carries full provenance in its HDF5 root
  attributes (scenario, scanner, crystal, budget, dose, master_seed, shard index, detector windows,
  n_phi/n_z). `config.toml` is the exact recipe: regenerate shard N with the base config +
  `--realization N`.

---

## The file contract

| file | what it is | who reads it |
|---|---|---|
| `lors_shardNNN.h5` | one shard's LOR list (schema `docs/SCHEMA.md`) + full provenance attrs | `thin_lm` → recon |
| `config.toml` | the run recipe (base config; shards = realizations 0..N) | regeneration; provenance |
| `scanner_geometry.json` | ring/panel geometry (+ crystal map if heterogeneous) | MLEM system model |
| `phantom/phantom_regions.csv` + `material_*.csv` | the phantom medium | build the μ-map for AC |
| `truth/depth_dose.csv` | dose(z), the SOBP distal edge | dose-R80 (the clinical range) |
| `truth/activity_profile_<budget>.csv` | binned true activity(z), per isotope + `total` | activity-R50 (the recon target) |
| `truth/{sobp_layers,run_meta,sampling_budget_<budget>}.csv` (+`_meta`) | beam design + normalization + per-isotope N_expected | scenario characterization |
| `README.md` | this contract | anyone browsing the tree |
| `SCHEMA.md` | the shard column schema (copy of `docs/SCHEMA.md`, refreshed each publish) | anyone reading the shards |
| `scanner_prods.pdf` | the productions note: the two scanners, constants, CTR calibration, statistics, column semantics | CryspBrainSim orientation |

`lors_shardNNN.h5` truth flag: `0` true, `1` scatter, `2` random. First-pass analysis is
**trues-only** (flag 0).

---

## Naming conventions

- **Shard files:** `lors_shardNNN.h5`, zero-padded (`shard000`…`shard009`) so they sort and glob
  cleanly and extend past 10 if ever needed.
- **Dose in the leaf:** `1Gy`, `0p5Gy` (decimal point → `p`).
- **Scanner name:** the `scanner.name` in `scanner_geometry.json`, crystal-neutral (e.g.
  `crysp_ring_1m` — NOT `crysp_csi_1m`; the crystal is a separate axis).
- **Crystal:** lower-case material key (`bgo`, `csi`).

---

## Engine gates (slots that exist but cannot be filled yet)

- `open_dualhead/` — needs partial-ring (`phi_gaps`) / planar-panel geometry in the simulator.
- mixed-crystal scanners (e.g. `bgo_core_csi_wings/`) — needs per-block crystal materials.
Closed-ring / head / children scanners are `CylShell` variants and are producible now. Expect only
producible arms in the tree.

---

*Authoritative copy: this file (with `publish_prod` in `PTCryspMC.jl`). `publish_prod` stamps a copy
at `PtCryspProds/README.md` so the tree is self-describing. The analysis repo references this
contract; it does not re-document it.*
