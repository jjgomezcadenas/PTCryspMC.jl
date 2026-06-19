module PTCryspMC

# Fast PET detector simulation: read a scenario (the annihilation source) and
# produce a coincidence list. See docs/pet_simulation.tex for the method and
# CLAUDE.md for the decisions and build order. dev/dev_steps.md logs the build.

using Random
using JSON

include("nist_data.jl")    # XCOM loader + interpolation
include("materials.jl")    # Material + macroscopic cross sections
include("geometry.jl")     # Cylinder + ray distance + phantom loader
include("sampling.jl")     # distance / process / Compton samplers + interaction kernel
include("transport.jl")    # photon-only transport through one volume
include("navigator.jl")    # multi-volume navigation across the geometry
include("source.jl")       # emission source + back-to-back acollinearity
include("singles.jl")      # chunked singles generation core (multi-thread driver + tests)
include("detector.jl")     # detector response: energy + position smearing
include("config.jl")       # TOML run-config reader (driver scripts)

export XCOMData, load_xcom,
       Material, load_material, load_materials, sigma_macro, mfp,
       Solid, Cylinder, Box, CylShell, Sphere, LogicalVolume, PhysicalVolume, Geometry,
       Scanner, r_outer, block_index, block_id, nblocks,
       solid, material, name, volume, mass,
       load_solid, load_geometry,
       is_inside, distance_to_exit, distance_to_entry,
       sample_interaction, sample_interaction_t, rotate_to_global_t,
       Interaction, Transported, propagate_photon,
       NavStep, locate, next_boundary, navigate_photon, navigate_single_photons,
       Source, UniformVolumeSource, PointSource, sample_point_in, sample_position, emit_pair, rand_direction,
       chunk_ranges, singles_chunk!,
       energy_fwhm, energy_sigma, smear_energy, smear_position,
       read_config, cfg_get, run_tag

end # module PTCryspMC
