# xsection-weighted-lors — status: deferred

The branch was created to let recorded coincidences carry their source-event
identity so cross-section replica variations could be applied by re-weighting
events (see ptcryspg4 docs/shared_plan.tex).

Decision (2026-07-26): **deferred indefinitely**. The production-level
cross-section systematic (u_xs = +-0.13 mm on the distal edge, from the
1000-replica fold in ptcryspg4) stands as the quoted systematic:
reconstruction is smooth near the edge, so couplings act at second order on
a 0.13 mm band (<0.03 mm). The error budget quotes sigma_R (statistics,
shard/thinning) and u_xs separately. If a reconstructed-level replica
distribution is ever wanted, the design remains: second parent id on
randoms, bank source mode reading ptcryspg4 source_bank.csv, weights
evaluator, resampler.

Active detector-side work instead: sigma_R reruns on the frozen data-driven
scenario `ptcrysp-scenarios/scenarios/uniform_headep_sobp_1e8_dd` (nominal
curves; existing unweighted pipeline, no code changes) — reference protocol,
CsI timing scan, one AFOV spot-check. Full record:
ptcryspg4 workshop/xsections_phases.md.
