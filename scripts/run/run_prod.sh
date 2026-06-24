#!/usr/bin/env zsh
# Production LOR chain for one or more run configs, run SEQUENTIALLY — each run is multi-threaded
# (it already uses the machine), so configs take turns rather than fight over the cores. Writes the
# production stack to prod/<tag>/{singles, lors_truth, randoms, lors_det}.h5.
#
# This is the PRODUCTION launcher. Its dev counterpart, scripts/run/run_matrix.sh, drives the
# full-stack dev chain (simulate_phantom -> build_coincidences -> output/, single-threaded, fanned
# out in parallel). This one drives the production chain instead:
#   simulate_source_mt -> build_true_coincidences_from_singles -> build_randoms_from_singles -> reco_lors
#
# HDF5 (so the runs carry the provenance attrs) + a PINNED nchunks, so the singles are
# bit-reproducible from (config + seed) — they can be regenerated exactly after pruning.
#
# Usage (from anywhere):
#   scripts/run/run_prod.sh                                    # every runs/*.toml
#   scripts/run/run_prod.sh sphere_water_csi sphere_water_bgo  # named configs
#   NEV=1000000 THREADS=8 scripts/run/run_prod.sh sphere_water_csi   # override scale/threads
#
# Env: NEV (events, default 1e7), NCHUNKS (chunks, default 144), THREADS (-t, default 16).
# Logs: /tmp/ptcprod_<tag>.log per config. Data: prod/<tag>/*.h5.

set -u
cd ${0:A:h:h:h}                       # repo root (in scripts/run/)

NEV=${NEV:-10000000}
NCHUNKS=${NCHUNKS:-144}
THREADS=${THREADS:-16}

# Collect configs: the args (name, name.toml or runs/name.toml all accepted), else all.
if (( $# )); then
  configs=()
  for a in "$@"; do configs+=(runs/${a:t:r}.toml); done
else
  configs=(runs/*.toml(N))            # (N) = no error if none match
fi
(( ${#configs} )) || { print "no configs found in runs/"; exit 1 }
for cfg in $configs; do
  [[ -f $cfg ]] || { print "missing config: $cfg"; exit 1 }
done

print "prod chain: ${#configs} config(s)  NEV=$NEV NCHUNKS=$NCHUNKS THREADS=$THREADS  -> ${(j: :)${configs:t:r}}"
start=$SECONDS

for cfg in $configs; do
  tag=${cfg:t:r}
  dir=prod/$tag
  print "=== $tag ==="
  if {
    julia -t $THREADS --project=. scripts/simulate_source_mt.jl --config $cfg --nevents $NEV --nchunks $NCHUNKS --format hdf5 &&
    julia -t $THREADS --project=. scripts/build_true_coincidences_from_singles.jl --config $cfg --singles $dir/singles.h5 &&
    julia --project=. scripts/build_randoms_from_singles.jl --config $cfg --singles $dir/singles.h5 &&
    julia --project=. scripts/reco_lors.jl --config $cfg
  } >/tmp/ptcprod_$tag.log 2>&1; then
    grep -hE "kept|acceptance" /tmp/ptcprod_$tag.log | tail -2
  else
    print "ERROR -> /tmp/ptcprod_$tag.log"
  fi
done

print "\nall ${#configs} finished in $((SECONDS - start)) s wall\n"
for cfg in $configs; do
  tag=${cfg:t:r}
  printf "%-24s " $tag
  grep -h "acceptance" /tmp/ptcprod_$tag.log 2>/dev/null | tail -1 || print "ERROR -> /tmp/ptcprod_$tag.log"
done
