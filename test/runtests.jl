using PTCryspMC
using Test

const DATA_DIR = joinpath(@__DIR__, "..", "data")

@testset "PTCryspMC" begin
    @testset "materials loading" begin
        # Single-material loader: reads only the named entry.
        w = load_material(DATA_DIR, "Water")
        @test w.name == "Water"
        @test w.density == 1.0
        @test !isempty(w.E)                       # has an XCOM grid

        # Water total cross section at 511 keV ~ 0.096 cm^-1 (mfp ~ 10.4 cm).
        Σtot = sum(sigma_macro(w, 0.511))
        @test isapprox(Σtot, 0.096, atol = 0.002)
        @test isapprox(mfp(w, 0.511), 1.0 / Σtot)

        # XCOM grid spans the cut (10 keV) up to 10 MeV, the range we care about.
        @test isapprox(first(w.E), 0.010)
        @test isapprox(last(w.E),  10.0)

        # Pair (last in the tuple) is below threshold at 511 keV but turns on at
        # high energy; Compton dominates throughout.
        C511, Ph511, P511 = sigma_macro(w, 0.511)
        @test P511 == 0.0
        @test C511 > Ph511                       # Compton >> photoelectric at 511 keV
        C10, Ph10, P10 = sigma_macro(w, 10.0)
        @test P10 > 0.0                           # pair contributes at 10 MeV
        @test C10 > P10 > Ph10                    # ordering of channels at 10 MeV

        # Vacuum: no XCOM table, zero cross sections, infinite mfp.
        v = load_material(DATA_DIR, "Vacuum")
        @test v.name == "Vacuum"
        @test isempty(v.E)
        @test sigma_macro(v, 0.511) == (0.0, 0.0, 0.0)
        @test mfp(v, 0.511) == Inf

        # Batch loader: builds every non-"_" entry, and agrees with the
        # single-material path for a shared material.
        ms = load_materials(DATA_DIR)
        @test haskey(ms, "Water") && haskey(ms, "Vacuum")
        @test !any(k -> startswith(k, "_"), keys(ms))   # doc keys skipped
        @test sigma_macro(ms["Water"], 0.511) == sigma_macro(w, 0.511)

        # An unknown name is an error, not a silent miss.
        @test_throws ErrorException load_material(DATA_DIR, "Nonexistent")
    end

    @testset "cylinder solid" begin
        c = Cylinder(8.0, 8.0)                    # R = H = 8 cm, centred at origin

        @test c isa Solid
        @test isapprox(volume(c), π * 8.0^2 * 16.0)   # πR²·(2H)

        @test is_inside(c, (0.0, 0.0, 0.0))
        @test is_inside(c, (8.0, 0.0, 8.0))       # on the rim: boundary is inside
        @test !is_inside(c, (8.01, 0.0, 0.0))

        # Interior pencil along +z exits through the far cap at 2H.
        @test isapprox(distance_to_exit((0.0, 0.0, -8.0), (0.0, 0.0, 1.0), c), 16.0)
        # Interior ray along +x exits through the lateral wall at R.
        @test isapprox(distance_to_exit((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), c), 8.0)
        # No entry from inside: the one forward crossing is an exit.
        @test distance_to_entry((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), c) == Inf

        # Exterior ray crossing the body: enters at the near wall, exits at the far.
        @test isapprox(distance_to_entry((-20.0, 0.0, 0.0), (1.0, 0.0, 0.0), c), 12.0)
        @test isapprox(distance_to_exit((-20.0, 0.0, 0.0), (1.0, 0.0, 0.0), c), 28.0)

        # A ray that misses the cylinder entirely: no entry, no exit.
        @test distance_to_entry((-20.0, 20.0, 0.0), (1.0, 0.0, 0.0), c) == Inf
        @test distance_to_exit((-20.0, 20.0, 0.0), (1.0, 0.0, 0.0), c) == Inf

        # Corner case: a ray aimed exactly at the rim (R, 0, H) must not slip
        # through the inclusive boundary and return Inf.
        @test isfinite(distance_to_exit((0.0, 0.0, 0.0), (8.0, 0.0, 8.0), c))
    end

    @testset "logical volume" begin
        w  = load_material(DATA_DIR, "Water")
        lv = LogicalVolume("phantom", Cylinder(8.0, 8.0), w)

        @test name(lv) == "phantom"
        @test solid(lv) === lv.solid
        @test material(lv) === w
        @test isapprox(volume(lv), π * 8.0^2 * 16.0)
        # Water Ø16×16 cm: mass = density · volume (≈ 3.2 kg).
        @test isapprox(mass(lv), 1.0 * π * 8.0^2 * 16.0)
    end

    @testset "physical volume placement" begin
        w   = load_material(DATA_DIR, "Water")
        lv  = LogicalVolume("phantom", Cylinder(8.0, 8.0), w)
        pv  = PhysicalVolume(lv, (10.0, 0.0, 0.0))   # placed off the origin

        # Accessors reach through the levels.
        @test solid(pv) === lv.solid
        @test material(pv) === w
        @test name(pv) == "phantom"
        @test isapprox(mass(pv), mass(lv))

        # The placement transform shifts the whole body by +10 in x.
        @test is_inside(pv, (10.0, 0.0, 0.0))         # centre, now at (10,0,0)
        @test !is_inside(pv, (0.0, 0.0, 0.0))         # old origin is now outside
        # Interior ray along +x exits the far wall at x = 18 → distance 8.
        @test isapprox(distance_to_exit((10.0, 0.0, 0.0), (1.0, 0.0, 0.0), pv), 8.0)
        # Exterior ray from the origin enters the near wall at x = 2 → distance 2.
        @test isapprox(distance_to_entry((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), pv), 2.0)
    end

    @testset "geometry loading" begin
        mats = load_materials(DATA_DIR)
        geom = load_geometry(joinpath(@__DIR__, "..", "geometry", "geometry.json"), mats)

        @test geom isa Geometry
        pv = geom.phantom
        @test pv isa PhysicalVolume
        @test name(pv) == "phantom"
        @test material(pv).name == "Water"
        @test isapprox(solid(pv).radius_cm, 8.0)
        @test isapprox(solid(pv).half_length_cm, 8.0)
        @test pv.position == (0.0, 0.0, 0.0)          # default = scanner centre

        # A non-origin placement is read from the phantom section's position_cm.
        tmp = tempname() * ".json"
        write(tmp, """{"phantom":{"shape":"cylinder","radius_cm":8.0,"half_length_cm":8.0,
                       "position_cm":[1.0,2.0,3.0],"material":"Water"}}""")
        @test load_geometry(tmp, mats).phantom.position == (1.0, 2.0, 3.0)
        rm(tmp)

        # A geometry file with no phantom section is an error.
        empty = tempname() * ".json"
        write(empty, """{"_doc":"nothing here"}""")
        @test_throws ErrorException load_geometry(empty, mats)
        rm(empty)

        # An unsupported shape is rejected, not silently loaded as a cylinder.
        bad = tempname() * ".json"
        write(bad, """{"phantom":{"shape":"sphere","radius_cm":5.0,"material":"Water"}}""")
        @test_throws ErrorException load_geometry(bad, mats)
        rm(bad)

        # A material the scene doesn't have is an error, not a silent miss.
        nomat = tempname() * ".json"
        write(nomat, """{"phantom":{"shape":"cylinder","radius_cm":8.0,"half_length_cm":8.0,"material":"Unobtainium"}}""")
        @test_throws ErrorException load_geometry(nomat, mats)
        rm(nomat)
    end
end
