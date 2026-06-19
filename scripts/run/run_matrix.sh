#!/usr/bin/env zsh
# Generate the LOR data (simulate_phantom -> build_coincidences) for one or more run
# configs, IN PARALLEL — one background job per config. Each config is independent (its
# own output/<tag>/ dir), so they fan out across the machine's cores. Plotting is a
# separate step: run scripts/run/plot_all.sh afterward (so you can re-plot without re-running
# the transport).
#
# Usage (from anywhere):
#   scripts/run/run_matrix.sh                                  # every runs/*.toml
#   scripts/run/run_matrix.sh sphere_water_csi cylinder_water_bgo   # named configs
#   scripts/run/run_matrix.sh runs/sphere_water_csi.toml            # paths also work
#
# Logs: /tmp/ptc_<tag>.log per config. Data: output/<tag>/{stack.csv, coincidences_*.csv}.

set -u
cd ${0:A:h:h:h}                       # repo root (in scripts/run/)

# Collect configs: the args (name, name.toml or runs/name.toml all accepted), else all.
if (( $# )); then
  configs=()
  for a in "$@"; do configs+=(runs/${a:t:r}.toml); done
else
  configs=(runs/*.toml(N))          # (N) = no error if none match
fi
(( ${#configs} )) || { print "no configs found in runs/"; exit 1 }

for cfg in $configs; do
  [[ -f $cfg ]] || { print "missing config: $cfg"; exit 1 }
done

print "launching ${#configs} config(s) in parallel: ${(j: :)${configs:t:r}}"
start=$SECONDS

for cfg in $configs; do
  tag=${cfg:t:r}
  (
    julia --project=. scripts/simulate_phantom.jl  --config $cfg &&
    julia --project=. scripts/build_coincidences.jl --config $cfg
  ) >/tmp/ptc_$tag.log 2>&1 &
done
wait

print "\nall ${#configs} finished in $((SECONDS - start)) s wall\n"
for cfg in $configs; do
  tag=${cfg:t:r}
  printf "%-24s " $tag
  grep -h "truth split" /tmp/ptc_$tag.log 2>/dev/null || print "ERROR -> /tmp/ptc_$tag.log"
done
