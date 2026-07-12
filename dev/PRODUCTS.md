# PtCryspProds — the products directory contract

The layout of the list-mode products produced by `PTCryspMC.jl` and consumed by the downstream
reconstruction / range-analysis repo. This is the interface between the two: the simulator writes
this tree (via `publish_prod`), the analysis repo reads it. For *why* the data is shaped this way
(shards, thinning, σ_R), read `dev/data_generation_strategy.md` first.

`PtCryspProds/` is a directory **sibling to the repos** (like `ptcrysp-scenarios/`), not committed
to `PTCryspMC.jl` (the `.h5` files are large and regenerable). It is the authoritative store of the
simulated LOR data.

---

## Generation (read this first)

Current products are **generation v2** — every shard carries `generation = "v2"` in its root
attributes, a **mandatory guard**: a consumer must check it and refuse to mix generations, because
the meaning of `t_decay_s` changed (see below). v2 vs the legacy off-centre masters:

- the phantom is **tumour-centred** (the SOBP dose-target centre at the ring centre — a fixed
  anatomical reference, `center_on="tumour"`, `source_z_offset_mm` stamped);
- timing is a family of **fixed acquisition scenarios** on the **irradiation-end clock** — the leaf
  is `del<t_del>s_ac<t_ac>s_<dose>` (arrival delay `t_del`, fixed length `t_ac`), **one leaf per
  scenario**, replacing the legacy `<budget>_<dose>`. `t_decay_s` has zero = **irradiation end**
  (attr `t_decay_zero`), not acquisition start;
- each LOR carries its **emitting isotope** (`isotope` column) — exact per-species work downstream;
- every shard is **self-describing**: the window (`t_del_s`/`t_ac_s`/`t1_s`/`t2_s`/`t_irr_s`), the
  full geometry embedded (`geometry_json`), and the Mizuno brain washout survival `washout_g` per
  isotope — **computed, NOT applied** (apply it downstream; see `washout_brain.tex`).

The authoritative v2 spec is **`dev/generation2_plan.md`**; this file is the directory layout.
Legacy off-centre `<budget>_<dose>` masters remain in the tree under their old token and never
collide with v2 leaves.

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
        del<t_del>s_ac<t_ac>s_<dose>/    ── ACQUISITION SCENARIO (v2): arrival delay t_del + fixed
                                            length t_ac + the master's TOP dose. One leaf per
                                            scenario. (Legacy off-centre: <budget>_<dose>.)
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
  crysp_ring_1m_csi_2x0/
    scanner_geometry.json
    csi/del120s_ac300s_1Gy/{config.toml, lors_shard000.h5 … lors_shard009.h5}
    csi/del180s_ac300s_1Gy/{config.toml, lors_shard000.h5 … }
    csi/del300s_ac300s_1Gy/{config.toml, lors_shard000.h5 … }
```
(v2: three acquisition-scenario leaves per crystal. Legacy off-centre example:
`crysp_ring_1m/{bgo,csi}/fast_1Gy/…`.)

---

## The axes

| level | axis | varies | shared with… |
|---|---|---|---|
| `<scenario>/` | proton field + phantom | scenario | the **source**, the **μ-map** (`phantom/`) & the **truth** (`truth/`) — across everything below |
| `<scanner>/` | scanner geometry | scanner | `scanner_geometry.json` — across the crystals in it |
| `<crystal>/` | crystal material (homogeneous only) | crystal | — |
| `del<t_del>s_ac<t_ac>s_<dose>/` | acquisition scenario (delay + length) + master top dose | scenario | — (v2; legacy: `<budget>_<dose>`) |
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
- **Acquisition scenario** (v2, the `del<t_del>s_ac<t_ac>s` leaf) = a fixed acquisition window
  `[t_del, t_del+t_ac]` on the irradiation-end clock; different scenarios are genuinely different
  sources (different N_j and randoms) — you cannot thin between them. (Legacy **budget** = fast /
  inroom / offline played the same role.) Within one scenario leaf, the 10 shards pool + thin as
  usual.
- **Dose** in the leaf name = the master's **top** dose. Lower doses are produced downstream by
  thinning *down*; they are not stored.

---

## Reading rules (for the analysis repo)

- **The master for one configuration:** `glob(<scenario>/<scanner>/<crystal>/del<t_del>s_ac<t_ac>s_<dose>/lors_shard*.h5)`
  → pool all shards → `thin_lm` (Bernoulli p = target/M_total over the union; see
  `dev/data_generation_strategy.md` §4) → realizations → σ_R. One σ_R curve **per acquisition
  scenario** (sweep the `del…` leaves to see the delay/washout dependence).
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
  attributes (scenario, scanner, crystal, dose, master_seed, shard index, detector windows,
  n_phi/n_z; **v2 adds** `generation`, `t_decay_zero`, `center_on`/`source_z_offset_mm`,
  `t_del_s`/`t_ac_s`/`t1_s`/`t2_s`/`t_irr_s`, `isotope_names`/`isotope_half_lives`, `geometry_json`,
  and the `washout_*`/`washout_g` set — the full enumeration is in `docs/SCHEMA.md`). `config.toml`
  is the exact recipe: regenerate shard N with the base config + `--realization N`.

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
**trues-only** (flag 0). v2 also carries the `isotope` column (0=O15, 1=C11, 2=N13, 3=C10, 4=O14;
gamma 1's decay for randoms) — for exact per-species selection and σ_R^(i); ignore it to emulate
the isotope-blind detector.

---

## Naming conventions

- **Shard files:** `lors_shardNNN.h5`, zero-padded (`shard000`…`shard009`) so they sort and glob
  cleanly and extend past 10 if ever needed.
- **Acquisition-scenario leaf (v2):** `del<t_del>s_ac<t_ac>s_<dose>` — integer seconds, e.g.
  `del120s_ac300s_1Gy`. (Legacy off-centre: `<budget>_<dose>`, e.g. `fast_1Gy`.)
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
