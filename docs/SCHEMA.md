# Output schema — PTCryspMC.jl

**Generated from the code by `scripts/gen_schema.jl` — do not edit by hand.** Regenerate with `julia --project=. scripts/gen_schema.jl`. Columns and types are introspected from `singles_columns` / `coinc_columns` and the struct field types; units and meaning come from the `singles_doc` / `coinc_doc` maps beside them. `test/runtests.jl` fails if this file drifts from the code.

Positions and energies are stored as quantized `Int16`: position = `round(mm / 0.1)`, energy = `round(keV / 0.1)` — lossless at detector resolution; decode with `decode_xyz` / `decode_e`. Times are `Float32` nanoseconds **relative to the decay**; the decay's absolute time in the acquisition window is the LOR column `t_decay_s` (`Float32` seconds, zero = acquisition start), so absolute photon time = `t_decay_s·1e9 + t`.

## Singles — `prod/<tag>/singles.{h5,csv}`

One row per detected photon (deposited in the ring; a miss writes nothing).

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `event_number` | Int32 |  | annihilation (decay) index |
| `gamma` | Int8 |  | which photon of the back-to-back pair (1 or 2) |
| `x_mm` | Int16 | mm | first crystal interaction point (the LOR point), x |
| `y_mm` | Int16 | mm | first crystal interaction point, y |
| `z_mm` | Int16 | mm | first crystal interaction point, z |
| `e_keV` | Int16 | keV | summed energy in the block (truth, unsmeared) |
| `iz` | Int16 |  | wheel (z) block index |
| `iphi` | Int16 |  | azimuthal (φ) block index |
| `nblocks` | Int8 |  | distinct blocks touched (1 = contained, >1 = overspill) |
| `n_scatter` | Int8 |  | phantom-scatter count (0 clean, 1 single, ≥2 multiple) |
| `isotope` | Int8 |  | emitter isotope id (API mode; 0 = single toy isotope otherwise) |
| `x0_mm` | Int16 | mm | annihilation (emission) point, x |
| `y0_mm` | Int16 | mm | annihilation point, y |
| `z0_mm` | Int16 | mm | annihilation point, z |
| `t_rel_ns` | Float32 | ns | photon time relative to its decay = TOF + scintillation jitter |

## LORs — `prod/<tag>/{lors_truth,lors_det,randoms}.h5`

One row per coincidence (line of response). All three files share this schema: `lors_truth` = same-annihilation pairs (truth, unsmeared); `lors_det` = the smeared, energy/DT-selected detector list (truth ∪ randoms); `randoms` = cross-decay accidentals (`truth = 2`, `dt_ns = NaN`).

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `event` | Int32 |  | annihilation index (shared by the two gammas; a cross-decay pair for randoms) |
| `truth` | Int8 |  | 0 = true, 1 = scatter, 2 = random |
| `x1_mm` | Int16 | mm | gamma 1 hit position (smeared in lors_det), x |
| `y1_mm` | Int16 | mm | gamma 1 hit position, y |
| `z1_mm` | Int16 | mm | gamma 1 hit position, z |
| `e1_keV` | Int16 | keV | gamma 1 energy (smeared in lors_det) |
| `t1_ns` | Float32 | ns | gamma 1 time relative to the decay |
| `iz1` | Int16 |  | gamma 1 wheel (z) block index |
| `iphi1` | Int16 |  | gamma 1 azimuthal (φ) block index |
| `nscat1` | Int8 |  | gamma 1 phantom-scatter count (0 clean, 1 single, ≥2 multiple) |
| `x2_mm` | Int16 | mm | gamma 2 hit position, x |
| `y2_mm` | Int16 | mm | gamma 2 hit position, y |
| `z2_mm` | Int16 | mm | gamma 2 hit position, z |
| `e2_keV` | Int16 | keV | gamma 2 energy |
| `t2_ns` | Float32 | ns | gamma 2 time relative to the decay |
| `iz2` | Int16 |  | gamma 2 wheel (z) block index |
| `iphi2` | Int16 |  | gamma 2 azimuthal (φ) block index |
| `nscat2` | Int8 |  | gamma 2 phantom-scatter count |
| `dt_ns` | Float32 | ns | timing residual (t1−t2) − TOF_diff; NaN for randoms |
| `x0_mm` | Int16 | mm | annihilation point (gamma 1's decay for randoms), x |
| `y0_mm` | Int16 | mm | annihilation point, y |
| `z0_mm` | Int16 | mm | annihilation point, z |
| `t_decay_s` | Float32 | s | absolute decay time within the acquisition window, zero = acquisition start (gamma 1's decay for randoms) |

### Truth code

| Value | Meaning |
|-------|---------|
| 0 | true |
| 1 | scatter |
| 2 | random |

A `truth = 1` (scatter) LOR is **single** scatter when `nscat1 + nscat2 == 1`, **multiple** when `≥ 2`.

## HDF5 root attributes (provenance)

Set by the writers (not part of the column lists), carried for provenance and exact regeneration of pruned singles. Representative keys:

| Attribute | Files | Meaning |
|-----------|-------|---------|
| `scenario_tag` | all | run tag (config filename) |
| `nrows` | all | row count |
| `xyz_scale_mm`, `e_scale_keV` | all | quantization scales (0.1 mm, 0.1 keV) |
| `seed` | singles | transport seed; on LORs, the builder's smear/time seed |
| `transport_seed`, `nchunks`, `nevents` | singles, lors_truth, lors_det | the recipe to regenerate the singles exactly |
| `crystal`, `phantom_material` | all | materials |
| `mode`, `has_randoms` | LORs | `truth`/`det`/`randoms`; whether randoms are merged |
| `tau_ns` | lors_det, randoms | coincidence window |
| `sigma_xyz_mm`, `eres`, `emin_keV`, `window_fwhm` | lors_det | detector response + energy cut |
| `t_relative_to_decay`, `t0_s`, `t1_s`, `half_life_s`, `time_seed` | LORs | the activity clock for absolute time |
| `t_decay_zero` | LORs | zero point of `t_decay_s` (`acquisition_start`) |
