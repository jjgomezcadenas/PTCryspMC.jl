module PTCryspMC

# Fast PET detector simulation: read a scenario (the annihilation source) and
# produce a coincidence list. See docs/pet_simulation.tex for the method and
# CLAUDE.md for the decisions and build order. dev/dev_steps.md logs the build.

using Random
using JSON

include("nist_data.jl")    # XCOM loader + interpolation
include("materials.jl")    # Material + macroscopic cross sections
include("geometry.jl")     # Cylinder + ray distance + phantom loader
include("sampling.jl")     # distance / process / Compton samplers
include("transport.jl")    # photon-only transport through a cylinder

export XCOMData, load_xcom,
       Material, load_material, load_materials, sigma_macro, mfp,
       Solid, Cylinder, Box, LogicalVolume, PhysicalVolume, Geometry,
       solid, material, name, volume, mass,
       load_solid, load_geometry,
       is_inside, distance_to_exit, distance_to_entry,
       Interaction, propagate_photon

end # module PTCryspMC
