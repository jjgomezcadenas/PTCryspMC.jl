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

export XCOMData, load_xcom,
       Material, load_material, load_materials, sigma_macro, mfp,
       Solid, Cylinder, Box, CylShell, Sphere, LogicalVolume, PhysicalVolume, Geometry,
       Scanner, r_outer, block_index, block_id, nblocks,
       solid, material, name, volume, mass,
       load_solid, load_geometry,
       is_inside, distance_to_exit, distance_to_entry,
       sample_interaction,
       Interaction, Transported, propagate_photon,
       NavStep, locate, next_boundary, navigate_photon,
       Source, UniformVolumeSource, PointSource, sample_point_in, sample_position, emit_pair, rand_direction

end # module PTCryspMC
