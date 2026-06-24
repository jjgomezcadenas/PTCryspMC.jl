using PTCryspMC
using Test
using Random

const DATA_DIR = joinpath(@__DIR__, "..", "data")
const GEOM_JSON = joinpath(@__DIR__, "..", "geometry", "geometry.json")

# Allocation probe for `navigate_single_photons`, kept at top level (a normal method, not a
# testset-local closure) and fed FIXED pos/dir so the measurement isolates the reducer. The rng
# advances each call, varying the histories. (`emit_pair` is also alloc-free now — its own probe is
# `_emit_alloc_run` below — but fixed pos/dir still keeps this measurement to just the reducer.)
function _singles_alloc_run(geom, rng, n)
    acc = 0.0
    pos = (0.0, 0.0, 0.0); dir = (1.0, 0.0, 0.0)   # phantom centre, radially into the ring
    for _ in 1:n
        acc += navigate_single_photons(geom, 0.511, pos, dir, rng).e
    end
    acc
end

# Allocation probe for `emit_pair` (the source path): the acollinearity tilt routes through the
# tuple twin `rotate_to_global_t`, so the whole annihilation draw is heap-free per call.
function _emit_alloc_run(src, rng, n)
    acc = 0.0
    for _ in 1:n
        _, d1, d2 = emit_pair(src, rng; acol_fwhm_deg=0.5)
        acc += d1[1] + d2[3]
    end
    acc
end

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

    @testset "CsI crystal material" begin
        csi = load_material(DATA_DIR, "CsI")
        @test csi.density == 4.51
        @test isapprox(first(csi.E), 0.010)
        @test isapprox(last(csi.E),  10.0)

        # Attenuation length at 511 keV ≈ 2.44 cm (CsI is a dense scintillator);
        # Compton dominates, photoelectric is well below it but non-negligible.
        C, Ph, P = sigma_macro(csi, 0.511)
        @test isapprox(1.0 / (C + Ph + P), 2.44, atol = 0.05)
        @test C > Ph > 0.0
        @test P == 0.0                                  # pair below threshold

        # The iodine K-edge (33.17 keV) must survive the duplicate-energy rows:
        # photoelectric jumps up sharply just above the edge vs. just below.
        Ph_below = sigma_macro(csi, 0.0330)[2]          # below the I K-edge
        Ph_above = sigma_macro(csi, 0.0340)[2]          # above it (below Cs edge)
        @test Ph_above > 2.0 * Ph_below
    end

    @testset "crystal scintillation properties (from the material DB)" begin
        # The crystal name pulls its intrinsic scintillation/readout properties — no config cards.
        csi = load_material(DATA_DIR, "CsI")            # cryogenic pure CsI
        @test csi.light_yield == 1.0e5
        @test csi.scint_decay_ns == [1000.0] && csi.scint_decay_w == [1.0]
        @test csi.eres_a == 0.05 && csi.pde == 0.45

        bgo = load_material(DATA_DIR, "BGO")            # cryogenic BGO — two decay components
        @test bgo.light_yield == 2.9e4
        @test bgo.scint_decay_ns == [747.0, 2253.0]
        @test bgo.scint_decay_w == [0.34, 0.66] && isapprox(sum(bgo.scint_decay_w), 1.0)
        @test bgo.eres_a == 0.10 && bgo.pde == 0.45

        # Non-scintillators get zero/empty defaults.
        w = load_material(DATA_DIR, "Water")
        @test w.light_yield == 0.0 && isempty(w.scint_decay_ns) && isempty(w.scint_decay_w) && w.eres_a == 0.0 && w.pde == 0.0

        # Loader validation: decay/weight length mismatch, and weights not summing to 1.
        @test_throws ErrorException PTCryspMC._build_material(DATA_DIR, "bad",
            Dict{String,Any}("density" => 1.0, "scint_decay_ns" => [1.0, 2.0], "scint_decay_w" => [1.0]))
        @test_throws ErrorException PTCryspMC._build_material(DATA_DIR, "bad",
            Dict{String,Any}("density" => 1.0, "scint_decay_ns" => [1.0, 2.0], "scint_decay_w" => [0.3, 0.3]))
    end

    @testset "BGO crystal material" begin
        bgo = load_material(DATA_DIR, "BGO")             # Bi4Ge3O12
        @test bgo.density == 7.13
        @test isapprox(first(bgo.E), 0.010)
        @test isapprox(last(bgo.E),  10.0)

        # Attenuation length at 511 keV ≈ 1.10 cm — much denser than CsI (2.44 cm);
        # high-Z Bi makes photoelectric large, though Compton still leads at 511 keV.
        C, Ph, P = sigma_macro(bgo, 0.511)
        @test isapprox(1.0 / (C + Ph + P), 1.10, atol = 0.05)
        @test C > Ph > 0.0
        @test P == 0.0                                  # pair below threshold

        # The bismuth K-edge (90.53 keV) must survive the duplicate-energy rows:
        # photoelectric jumps up sharply just above the edge vs. just below.
        Ph_below = sigma_macro(bgo, 0.0900)[2]          # below the Bi K-edge
        Ph_above = sigma_macro(bgo, 0.0910)[2]          # above it
        @test Ph_above > 2.0 * Ph_below
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

    @testset "box solid" begin
        b = Box(2.4, 2.4, 1.85)                   # the CsI crystal: 48 x 48 x 37 mm
        @test b isa Solid
        @test isapprox(volume(b), 8.0 * 2.4 * 2.4 * 1.85)

        @test is_inside(b, (0.0, 0.0, 0.0))
        @test is_inside(b, (2.4, -2.4, 1.85))     # corner: boundary is inside
        @test !is_inside(b, (2.41, 0.0, 0.0))

        # Interior ray along +z exits the far face at half_z = 1.85.
        @test isapprox(distance_to_exit((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), b), 1.85)
        # Interior ray along +x exits the side at half_x = 2.4.
        @test isapprox(distance_to_exit((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), b), 2.4)
        @test distance_to_entry((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), b) == Inf  # inside

        # Exterior ray entering the -z face and leaving the +z face.
        @test isapprox(distance_to_entry((0.0, 0.0, -5.0), (0.0, 0.0, 1.0), b), 3.15)
        @test isapprox(distance_to_exit((0.0, 0.0, -5.0), (0.0, 0.0, 1.0), b), 6.85)

        # A ray parallel to the box but outside it: no hit.
        @test distance_to_entry((10.0, 0.0, -5.0), (0.0, 0.0, 1.0), b) == Inf
        @test distance_to_exit((10.0, 0.0, -5.0), (0.0, 0.0, 1.0), b) == Inf

        # The shoot setup: box placed so its entry face is at world z = 0
        # (crystal spans z in [0, 3.7] cm). A photon entering there along +z
        # traverses the full 3.7 cm depth.
        pv = PhysicalVolume(LogicalVolume("crystal", b, load_material(DATA_DIR, "CsI")),
                            (0.0, 0.0, 1.85))
        @test is_inside(pv, (0.0, 0.0, 0.0))      # entry face
        @test is_inside(pv, (0.0, 0.0, 3.7))      # back face
        @test !is_inside(pv, (0.0, 0.0, 3.8))
        @test isapprox(distance_to_exit((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), pv), 3.7)
    end

    @testset "cylindrical shell solid" begin
        cs = CylShell(38.7, 3.7, 51.2)            # CRYSP1M ring: Ri=38.7, Ro=42.4
        @test cs isa Solid
        @test isapprox(r_outer(cs), 42.4)
        @test isapprox(volume(cs), π * (42.4^2 - 38.7^2) * 2 * 51.2)

        @test !is_inside(cs, (0.0, 0.0, 0.0))     # the bore
        @test is_inside(cs, (40.0, 0.0, 0.0))     # in the wall
        @test is_inside(cs, (38.7, 0.0, 0.0))     # inner wall (inclusive)
        @test is_inside(cs, (42.4, 0.0, 0.0))     # outer wall (inclusive)
        @test !is_inside(cs, (50.0, 0.0, 0.0))    # beyond the outer wall
        @test !is_inside(cs, (40.0, 0.0, 52.0))   # beyond the cap

        # From the bore centre along +x: first wall entry at the inner radius.
        @test isapprox(distance_to_entry((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), cs), 38.7)
        # From the inner wall (x=38.7) along +x: exit the outer wall 3.7 cm on.
        @test isapprox(distance_to_exit((38.7, 0.0, 0.0), (1.0, 0.0, 0.0), cs), 3.7)
        # In the bore, distance_to_exit is Inf (the point is not in the wall).
        @test distance_to_exit((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), cs) == Inf

        # The diametral chord from x=-50 along +x crosses outer, bore, bore, outer:
        # enters the near wall at x=-42.4 …
        @test isapprox(distance_to_entry((-50.0, 0.0, 0.0), (1.0, 0.0, 0.0), cs), 7.6)
        # … from inside that near wall (x=-40) it exits into the bore at x=-38.7 …
        @test isapprox(distance_to_exit((-40.0, 0.0, 0.0), (1.0, 0.0, 0.0), cs), 1.3)
        # … and from the bore (x=-20) the next wall entry is the far inner wall x=38.7.
        @test isapprox(distance_to_entry((-20.0, 0.0, 0.0), (1.0, 0.0, 0.0), cs), 58.7)

        # An axial ray inside the wall exits through a cap.
        @test isapprox(distance_to_exit((40.0, 0.0, 0.0), (0.0, 0.0, 1.0), cs), 51.2)
        # A ray straight down the bore axis never touches the wall.
        @test distance_to_entry((0.0, 0.0, -100.0), (0.0, 0.0, 1.0), cs) == Inf
    end

    @testset "sphere solid" begin
        s = Sphere(5.0)
        @test s isa Solid
        @test isapprox(volume(s), (4 / 3) * π * 5.0^3)

        @test is_inside(s, (0.0, 0.0, 0.0))       # centre
        @test is_inside(s, (5.0, 0.0, 0.0))       # surface (inclusive)
        @test is_inside(s, (3.0, 4.0, 0.0))       # on the surface (r=5)
        @test !is_inside(s, (5.1, 0.0, 0.0))      # just outside

        # From the centre, the exit is one radius away in any direction; no entry.
        @test isapprox(distance_to_exit((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), s), 5.0)
        @test isapprox(distance_to_exit((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), s), 5.0)
        @test distance_to_entry((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), s) == Inf

        # From outside heading at the centre: enter the near face, exit the far one.
        @test isapprox(distance_to_entry((-10.0, 0.0, 0.0), (1.0, 0.0, 0.0), s), 5.0)
        @test isapprox(distance_to_exit((-10.0, 0.0, 0.0), (1.0, 0.0, 0.0), s), 15.0)

        # Heading away, a miss (offset > R), and a tangent all give Inf.
        @test distance_to_entry((10.0, 0.0, 0.0), (1.0, 0.0, 0.0), s) == Inf
        @test distance_to_entry((-10.0, 6.0, 0.0), (1.0, 0.0, 0.0), s) == Inf
        @test distance_to_exit((10.0, 0.0, 0.0), (1.0, 0.0, 0.0), s) == Inf

        # Loaded from JSON like the other shapes.
        sl = load_solid(Dict("shape" => "sphere", "radius_cm" => 8.0))
        @test sl isa Sphere && isapprox(sl.radius_cm, 8.0)
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
        geom = load_geometry(GEOM_JSON, mats)
        @test geom isa Geometry

        # world: the Air mother volume, non-interacting, enclosing the daughters.
        @test name(geom.world) == "world"
        @test material(geom.world).name == "Air"
        @test sigma_macro(material(geom.world), 0.511) == (0.0, 0.0, 0.0)
        @test isapprox(solid(geom.world).radius_cm, 60.0)
        @test isapprox(solid(geom.world).half_length_cm, 60.0)

        # phantom: a daughter at the origin.
        pv = geom.phantom
        @test pv isa PhysicalVolume
        @test name(pv) == "phantom"
        @test material(pv).name == "Water"
        @test isapprox(solid(pv).radius_cm, 8.0)
        @test pv.position == (0.0, 0.0, 0.0)

        # scanner: the detector ring daughter.
        @test geom.scanner isa Scanner

        # build a geometry file with a valid world + the given phantom JSON object
        tmpgeo(phantom) = (p = tempname() * ".json";
            write(p, """{"world":{"shape":"cylinder","radius_cm":60,"half_length_cm":60,"material":"Air"},
                         "phantom":$phantom}"""); p)

        # A non-origin placement is read from the phantom section's position_cm.
        tmp = tmpgeo("""{"shape":"cylinder","radius_cm":8.0,"half_length_cm":8.0,"position_cm":[1.0,2.0,3.0],"material":"Water"}""")
        @test load_geometry(tmp, mats).phantom.position == (1.0, 2.0, 3.0)
        # the scanner section is optional: absent → nothing.
        @test load_geometry(tmp, mats).scanner === nothing
        rm(tmp)

        # Missing world section is an error.
        nw = tempname() * ".json"
        write(nw, """{"phantom":{"shape":"cylinder","radius_cm":8.0,"half_length_cm":8.0,"material":"Water"}}""")
        @test_throws ErrorException load_geometry(nw, mats)
        rm(nw)

        # Missing phantom section is an error.
        np = tempname() * ".json"
        write(np, """{"world":{"shape":"cylinder","radius_cm":60,"half_length_cm":60,"material":"Air"}}""")
        @test_throws ErrorException load_geometry(np, mats)
        rm(np)

        # An unsupported shape is rejected, not silently loaded as a cylinder.
        bad = tmpgeo("""{"shape":"torus","radius_cm":5.0,"material":"Water"}""")
        @test_throws ErrorException load_geometry(bad, mats)
        rm(bad)

        # A material the scene doesn't have is an error, not a silent miss.
        nomat = tmpgeo("""{"shape":"cylinder","radius_cm":8.0,"half_length_cm":8.0,"material":"Unobtainium"}""")
        @test_throws ErrorException load_geometry(nomat, mats)
        rm(nomat)
    end

    @testset "air world transport" begin
        mats = load_materials(DATA_DIR)
        air  = load_geometry(GEOM_JSON, mats).world      # the non-interacting mother
        rng  = MersenneTwister(1)

        # A photon in pure Air takes one straight step to the world boundary,
        # unchanged in energy, depositing nothing.
        recs = propagate_photon(0.511, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), air, rng).recs
        @test length(recs) == 1
        @test recs[1].process == :escape
        @test isapprox(recs[1].z, 60.0)        # exits the +z face at half_length
        @test isapprox(recs[1].e_in, 0.511)    # no energy loss
        @test recs[1].e_dep == 0.0

        # Along +x it exits the lateral wall at the world radius.
        recs2 = propagate_photon(0.511, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), air, rng).recs
        @test length(recs2) == 1 && recs2[1].process == :escape
        @test isapprox(recs2[1].x, 60.0)
        @test isapprox(recs2[1].e_in, 0.511)
    end

    @testset "scanner block / wheel grid" begin
        mats = load_materials(DATA_DIR)
        sc   = load_geometry(GEOM_JSON, mats).scanner

        @test sc isa Scanner
        @test sc.n_phi == 48 && sc.n_z == 20
        @test nblocks(sc) == 960
        @test name(sc.volume) == "CRYSP_CSI_1M"
        @test material(sc.volume).name == "CsI"
        @test solid(sc.volume) isa CylShell
        @test isapprox(solid(sc.volume).r_inner_cm, 38.7)

        # +x axis at the axial centre (z=0) → sector 0, the middle wheel 10 of 20.
        @test block_index(sc, (40.0, 0.0, 0.0)) == (0, 10)
        @test block_id(sc, (40.0, 0.0, 0.0)) == 10 * 48 + 0
        # +x axis near the −z end → block (0, 0), id 0.
        @test block_index(sc, (40.0, 0.0, -51.0)) == (0, 0)
        @test block_id(sc, (40.0, 0.0, -51.0)) == 0
        # φ = π (the −x side) → azimuthal sector 24 of 48.
        @test block_index(sc, (-40.0, 0.0, 0.0))[1] == 24
        # z near +H → wheel 19; the top edge clamps to 19.
        @test block_index(sc, (40.0, 0.0, 51.0))[2] == 19
        @test block_index(sc, (40.0, 0.0, 51.2))[2] == 19
        # linear id follows iz·n_phi + iφ.
        iφ, iz = block_index(sc, (-40.0, 0.0, 51.0))
        @test block_id(sc, (-40.0, 0.0, 51.0)) == iz * 48 + iφ
    end

    @testset "locate" begin
        mats = load_materials(DATA_DIR)
        geom = load_geometry(GEOM_JSON, mats)

        @test locate(geom, (0.0, 0.0, 0.0))   === geom.phantom          # phantom centre
        @test locate(geom, (5.0, 0.0, 0.0))   === geom.phantom          # inside phantom
        @test locate(geom, (20.0, 0.0, 0.0))  === geom.world            # air bore
        @test locate(geom, (40.0, 0.0, 0.0))  === geom.scanner.volume   # in the ring wall
        @test locate(geom, (50.0, 0.0, 0.0))  === geom.world            # air, beyond the ring
        @test locate(geom, (100.0, 0.0, 0.0)) === nothing               # escaped the world
    end

    @testset "next_boundary + bore re-entry" begin
        mats = load_materials(DATA_DIR)
        geom = load_geometry(GEOM_JSON, mats)

        # From the air bore heading +x at the ring: the near (inner) wall at 38.7.
        d = next_boundary(geom, geom.world, (20.0, 0.0, 0.0), (1.0, 0.0, 0.0))
        @test isapprox(d, 18.7)                          # 38.7 − 20

        # Bore re-entry capability: a chord that misses the central phantom (closest
        # approach 20 cm > the 8 cm phantom) crosses the bore and lands on the FAR
        # inner wall — the opposite crystal a backscattered photon would re-enter.
        start, dir = (-30.0, 20.0, 0.0), (1.0, 0.0, 0.0)
        d = next_boundary(geom, geom.world, start, dir)
        @test isfinite(d)
        far = (start[1] + d*dir[1], start[2] + d*dir[2], start[3])
        @test isapprox(hypot(far[1], far[2]), 38.7; atol=1e-6)   # the inner wall
        @test is_inside(geom.scanner.volume, far)
        @test far[1] > 0                                          # the OPPOSITE side
    end

    @testset "navigate_photon — reduction to single-volume" begin
        mats = load_materials(DATA_DIR)
        geom = load_geometry(GEOM_JSON, mats)

        # A photon launched in the air bore (outside the phantom) straight at the ring
        # must give the SAME crystal deposits as a bare single-volume transport from the
        # inner wall: the air leg consumes no randomness, so with an identical seed the
        # ring segment is bit-for-bit the same. This is the navigator reducing to
        # `propagate_photon` when one absorbing material is crossed.
        start, dir = (20.0, 0.0, 0.0), (1.0, 0.0, 0.0)   # r=20 > phantom (8): misses it
        steps = navigate_photon(geom, 0.511, start, dir, MersenneTwister(7))
        scan  = filter(s -> s.volume == :scanner, steps)
        @test !isempty(scan)

        entry = (38.7, 0.0, 0.0)                          # where the air skip lands
        ref = propagate_photon(0.511, entry, dir, geom.scanner.volume, MersenneTwister(7)).recs
        @test length(scan) >= length(ref)                 # ≥ (bore re-entry may add more)
        for i in eachindex(ref)                           # the first ring segment matches
            @test scan[i].hit.process == ref[i].process
            @test isapprox(scan[i].hit.e_dep, ref[i].e_dep)
            @test isapprox(scan[i].hit.x, ref[i].x; atol=1e-9)
            @test isapprox(scan[i].hit.y, ref[i].y; atol=1e-9)
        end
        # off the scanner the block tags are -1; on it they are valid.
        for s in steps
            if s.volume == :scanner
                @test 0 <= s.iphi < geom.scanner.n_phi && 0 <= s.iz < geom.scanner.n_z
            else
                @test s.iz == -1 && s.iphi == -1
            end
        end
    end

    @testset "navigate_photon — phantom leg + conservation" begin
        mats = load_materials(DATA_DIR)
        geom = load_geometry(GEOM_JSON, mats)
        rng  = MersenneTwister(11)

        # Monte Carlo property check over 300 photons born in the phantom: accumulate
        # the invariants and assert each ONCE, so the test count reflects the properties
        # checked, not the number of samples.
        saw_phantom = false
        conserved = true; phantom_inside = true; scanner_inside = true
        for _ in 1:300
            c = 2rand(rng) - 1; φ = 2π*rand(rng); s = sqrt(1 - c^2)
            dir = (s*cos(φ), s*sin(φ), c)
            steps = navigate_photon(geom, 0.511, (0.0, 0.0, 0.0), dir, rng)  # born in phantom

            conserved &= sum(st.hit.e_dep for st in steps) <= 0.511 + 1e-9   # across all volumes
            for st in steps
                if st.volume == :phantom && st.hit.process != :escape
                    saw_phantom = true
                    # phantom interactions lie inside the water cylinder (r≤8, |z|≤8).
                    phantom_inside &= hypot(st.hit.x, st.hit.y) <= 8.0 + 1e-6 &&
                                      abs(st.hit.z) <= 8.0 + 1e-6
                elseif st.volume == :scanner
                    # scanner interactions lie inside the ring wall (r in [38.7, 42.4]).
                    scanner_inside &= 38.7 - 1e-6 <= hypot(st.hit.x, st.hit.y) <= 42.4 + 1e-6
                end
            end
        end
        @test conserved          # energy conserved across volumes, every event
        @test phantom_inside     # every phantom interaction inside the water cylinder
        @test scanner_inside     # every scanner interaction inside the ring wall
        @test saw_phantom        # some photons Compton-scatter in the water phantom

        # A straight axial line misses the ring geometrically (r = 0 never reaches the
        # 38.7 cm inner wall).
        @test !isfinite(distance_to_entry((0.0, 0.0, 0.0), (0.0, 0.0, 1.0),
                                           geom.scanner.volume))

        # But a photon fired down the axis can still light the ring: it Compton-scatters
        # in the phantom, deflects off-axis, and the scattered photon reaches a crystal —
        # the patient-scatter path the navigator must follow across volumes. (Seed chosen
        # to scatter; energy is fully conserved phantom + ring = 511 keV.)
        steps = navigate_photon(geom, 0.511, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0),
                                MersenneTwister(2))
        @test any(s -> s.volume == :phantom && s.hit.process == :compton, steps)
        @test any(s -> s.volume == :scanner, steps)
        @test isapprox(sum(s.hit.e_dep for s in steps), 0.511; atol=1e-9)
    end

    @testset "rotate_to_global_t == rotate_to_global" begin
        # The tuple twin must reproduce the Vector version bit-for-bit (it is the basis of
        # the alloc-free singles path). Check across general and near-pole reference axes.
        rng = MersenneTwister(3)
        for _ in 1:500
            lv = (randn(rng), randn(rng), randn(rng))
            ax = (randn(rng), randn(rng), randn(rng))
            v = PTCryspMC.rotate_to_global(Float64[lv[1], lv[2], lv[3]], ax)
            t = rotate_to_global_t(lv[1], lv[2], lv[3], ax)
            @test t[1] == v[1] && t[2] == v[2] && t[3] == v[3]
        end
        for ax in ((0.0, 0.0, 1.0), (0.0, 0.0, -1.0), (1e-7, 0.0, 1.0))   # poles
            lv = (0.3, -0.5, 0.8)
            v = PTCryspMC.rotate_to_global(Float64[lv...], ax)
            t = rotate_to_global_t(lv..., ax)
            @test t[1] == v[1] && t[2] == v[2] && t[3] == v[3]
        end
    end

    @testset "navigate_single_photons — matches the full-stack reduction" begin
        mats = load_materials(DATA_DIR)
        geom = load_geometry(GEOM_JSON, mats)

        # navigate_single_photons folds a photon to the singles summary with no allocation;
        # it must return EXACTLY _reduce_steps(navigate_photon(...)) on the same history —
        # same physics, same RNG draw order. This guards the two shells against drift.
        # Both photons of back-to-back pairs from the phantom centre, so the phantom leg,
        # the ring, overspill and absorption are all exercised.
        nmismatch = 0; ntested = 0
        for ev in 1:3000
            r_emit = MersenneTwister(2000 + ev)
            pos, d1, d2 = emit_pair(PointSource(geom.phantom.position), r_emit)
            for dir in (d1, d2)
                rf = MersenneTwister(7000 + ev)              # identical streams for both calls
                rs = MersenneTwister(7000 + ev)
                full = PTCryspMC._reduce_steps(navigate_photon(geom, 0.511, pos, dir, rf))
                sing = navigate_single_photons(geom, 0.511, pos, dir, rs)
                ntested += 1
                full == sing || (nmismatch += 1)
            end
        end
        @test ntested == 6000
        @test nmismatch == 0

        # The reducer reaches the ring (the summary is non-trivial) and respects the
        # contract: a reached photon has a valid block, an unreached one is the empty single.
        r = MersenneTwister(42)
        reached_any = false; valid = true
        for ev in 1:2000
            pos, d1, d2 = emit_pair(PointSource(geom.phantom.position), r)
            s = navigate_single_photons(geom, 0.511, pos, d1, r)
            if s.reached
                reached_any = true
                valid &= (0 <= s.iphi < geom.scanner.n_phi) && (0 <= s.iz < geom.scanner.n_z) &&
                         s.nblocks >= 1 && s.e > 0.0
            else
                # Unreached → empty LOR fields (but phscat may be true: it can scatter in
                # the phantom and still never reach a crystal).
                valid &= s.x == 0.0 && s.y == 0.0 && s.z == 0.0 && s.iz == -1 &&
                         s.iphi == -1 && s.e == 0.0 && s.nblocks == 0
            end
        end
        @test reached_any
        @test valid
    end

    @testset "navigate_single_photons — zero allocation" begin
        mats = load_materials(DATA_DIR)
        geom = load_geometry(GEOM_JSON, mats)

        # The whole reason the function exists: fold a photon with no heap traffic, so GC
        # cannot serialise the threads at 10⁸ decays. Measure with `_singles_alloc_run`
        # (fixed inputs, only the reducer in the loop): per-photon allocation = the slope,
        # which must be exactly zero.
        _singles_alloc_run(geom, MersenneTwister(5), 10)     # compile
        b1000 = @allocated _singles_alloc_run(geom, MersenneTwister(5), 1000)
        b2000 = @allocated _singles_alloc_run(geom, MersenneTwister(5), 2000)
        @test (b2000 - b1000) == 0
    end

    @testset "simulate_source_mt — chunking & determinism" begin
        mats = load_materials(DATA_DIR)
        geom = load_geometry(GEOM_JSON, mats)
        src  = UniformVolumeSource(geom.phantom)

        # 1. Partition cover: the chunk ranges exactly tile 1:N, no gap/overlap — including
        #    nchunks > N (clamped) and nchunks = 1.
        for (N, k) in ((100, 8), (7, 16), (1000, 1), (5, 5), (333, 7))
            r = chunk_ranges(N, k)
            @test reduce(vcat, collect.(r)) == collect(1:N)   # exact tiling, in order
            @test length(r) == max(1, min(k, N))
        end

        # In-memory generation through the SAME core the driver uses, via a vector sink;
        # `threaded` toggles the @threads loop vs a sequential one over identical chunks/RNGs.
        function gen(; nevents, seed, nchunks, threaded)
            ranges = chunk_ranges(nevents, nchunks)
            rngs   = [MersenneTwister(seed + (c - 1)) for c in eachindex(ranges)]
            parts  = [Tuple[] for _ in ranges]
            body(c) = singles_chunk!(geom, src, 0.511, 0.010, 0.5, ranges[c], rngs[c], mats["CsI"]) do ev, g, s, _, t_rel
                push!(parts[c], (ev, g, s.iz, s.iphi, round(s.e, digits=9), s.nblocks, s.phscat, round(t_rel, digits=6)))
            end
            if threaded
                Threads.@threads for c in eachindex(ranges); body(c); end
            else
                for c in eachindex(ranges); body(c); end
            end
            reduce(vcat, parts)                               # glue in chunk order
        end

        # 2. Thread-invariance: the @threads loop produces EXACTLY the sequential result —
        #    no race, no shared RNG across chunks (the core MT-correctness guarantee).
        @test gen(nevents=4000, seed=1, nchunks=16, threaded=true) ==
              gen(nevents=4000, seed=1, nchunks=16, threaded=false)

        # 3. Event-ordered output (so the streaming build_coincidences reader needs no change).
        rows = gen(nevents=4000, seed=1, nchunks=16, threaded=true)
        @test !isempty(rows)
        @test issorted(rows; by=first)
    end

    @testset "singles quantization + binary part + HDF5 round-trip" begin
        # Quantization is lossless to ≤ half a step, and bounds-checked (fails loud, no wrap).
        for mm in (0.0, 12.34, -250.0, 511.9)
            @test abs(decode_xyz(encode_xyz_mm(mm)) - mm) <= XYZ_SCALE_MM/2 + 1e-9
        end
        for keV in (0.0, 1.0, 250.5, 511.0)
            @test abs(decode_e(encode_e_keV(keV)) - keV) <= E_SCALE_KEV/2 + 1e-9
        end
        @test_throws ErrorException encode_xyz_mm(4000.0)   # 40000 > Int16 @0.1mm
        @test_throws ErrorException encode_e_keV(-1.0)

        # A small synthetic buffer (positions cm, energy MeV — as navigate_single_photons returns).
        buf = SinglesBuffer()
        for i in 1:200
            s = (reached=true, x=0.13i, y=-0.07i, z=0.05i, iz=i % 20, iphi=i % 48,
                 e=0.001*(50 + i), nblocks=(i % 3) + 1, phscat=isodd(i))
            push_single!(buf, i, (i % 2) + 1, s, (0.01i, -0.02i, 0.03i), 0.05i)
        end
        @test length(buf) == 200

        # Binary part round-trip (the thread-safe intermediate).
        io = IOBuffer(); write_part(io, buf); seekstart(io)
        rb = read_part(io, 200)
        for ((_, va), (_, vb)) in zip(singles_columns(buf), singles_columns(rb))
            @test va == vb
        end

        # part-file → HDF5 pack → stream-read round-trip (multiple slices via small batch).
        dir = mktempdir()
        part = joinpath(dir, "p1.bin"); open(o -> write_part(o, buf), part, "w")
        h5 = joinpath(dir, "s.h5")
        @test pack_singles_hdf5(h5, [part], [200], Dict("seed" => 7)) == 200
        @test !isfile(part)                              # part consumed by the packer
        got = SinglesBuffer()
        foreach_singles_hdf5(h5; batch=64) do b
            for ((_, dst), (_, srcv)) in zip(singles_columns(got), singles_columns(b))
                append!(dst, srcv)
            end
        end
        @test length(got) == 200
        for ((_, va), (_, vb)) in zip(singles_columns(buf), singles_columns(got))
            @test va == vb
        end
        rm(dir; recursive=true)
    end

    @testset "coincidence core — singles fill ≡ full-stack fill" begin
        mats = load_materials(DATA_DIR)
        geom = load_geometry(GEOM_JSON, mats)

        # Reduce a full-stack NavStep history to (GammaAcc, phscat) the build_coincidences way.
        # `GammaAcc` holds only the accumulated geometry/energy; the phantom-scatter flag is carried
        # alongside (set by the caller, not the accumulator). Positions cm→mm (×10), energy MeV→keV.
        function full_acc(steps)
            a = GammaAcc(); ps = false
            for st in steps
                st.hit.process === :escape && continue
                if st.volume === :scanner
                    fill_full!(a, st.hit.x*10, st.hit.y*10, st.hit.z*10, st.hit.e_dep*1000, st.iz, st.iphi)
                elseif st.volume === :phantom
                    ps = true
                end
            end
            (a, ps)
        end
        sing_acc(s) = (a = GammaAcc(); s.reached && fill_singles!(a, s.x*10, s.y*10, s.z*10, s.e*1000, s.iz, s.iphi, s.nblocks); (a, s.phscat))
        # Compare the (acc, phscat) pairs for REACHED gammas (the only ones that can form a LOR):
        # discrete fields + LOR point EXACT, phantom-scatter flag equal; the summed energy agrees only
        # to float precision (full sums per-interaction keV, singles sums MeV once then ×1000).
        function feq(A, B)
            (af, aps), (bf, bps) = A, B
            af.reached == bf.reached && (!af.reached || (
                af.x == bf.x && af.y == bf.y && af.z == bf.z && af.iz == bf.iz && af.iphi == bf.iphi &&
                af.overspill == bf.overspill && aps == bps && isapprox(af.e, bf.e; rtol=1e-9)))
        end
        # A LOR tuple (20 fields: …,iz2(15),dt(16),x0(17),y0,z0,truth(20)). Discrete fields, the LOR
        # point, AND the timestamps/dt (6,13,16) match exactly — `t` is computed once and passed
        # identically to both fills' finish_event!. Only the summed energies (5,12) differ, at float
        # precision (full sums per-interaction keV, singles sums MeV once then ×1000).
        ceq(p, q) = all(j -> p[j] == q[j], (1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20)) &&
                    isapprox(p[5], q[5]; rtol=1e-9) && isapprox(p[12], q[12]; rtol=1e-9)

        # The singles fill (from navigate_single_photons) must produce the SAME GammaAcc as the
        # full-stack fill (from navigate_photon) — discrete fields exact, energy to float precision —
        # and hence the SAME coincidences through the shared finish_event!. The per-gamma time `t`
        # (TOF+jitter) and phantom-scatter flag are passed to finish_event! (not on the accumulator);
        # one `t` is computed and shared by both paths. TRUTH mode: no smearing, no energy cut.
        resp = Response(0.0, 0.0, 0.0, false, 0.0)
        mat  = mats["CsI"]; trng = MersenneTwister(1)
        emit_full = Any[]; emit_sing = Any[]
        nacc_mismatch = 0
        for ev in 1:3000
            pos, d1, d2 = emit_pair(PointSource(geom.phantom.position), MersenneTwister(5000 + ev))
            af = GammaAcc[]; as = GammaAcc[]; psf = Bool[]; pss = Bool[]; gt = Float64[]
            for (k, dir) in ((1, d1), (2, d2))
                seed = 90000 + 2ev + k
                gf2 = full_acc(navigate_photon(geom, 0.511, pos, dir, MersenneTwister(seed)))
                gs2 = sing_acc(navigate_single_photons(geom, 0.511, pos, dir, MersenneTwister(seed)))
                feq(gf2, gs2) || (nacc_mismatch += 1)
                t = gs2[1].reached ?                 # one per-gamma time, shared by both paths
                    tof_ns((pos[1]*10, pos[2]*10, pos[3]*10), (gs2[1].x, gs2[1].y, gs2[1].z)) +
                        first_photon_jitter(mat, gs2[1].e * 1e-3, trng) : 0.0
                push!(af, gf2[1]); push!(psf, gf2[2])
                push!(as, gs2[1]); push!(pss, gs2[2]); push!(gt, t)
            end
            finish_event!((a...) -> push!(emit_full, a), ev, af[1], af[2], gt[1], gt[2], psf[1], psf[2], pos[1], pos[2], pos[3], resp, MersenneTwister(ev))
            finish_event!((a...) -> push!(emit_sing, a), ev, as[1], as[2], gt[1], gt[2], pss[1], pss[2], pos[1], pos[2], pos[3], resp, MersenneTwister(ev))
        end
        @test nacc_mismatch == 0             # GammaAcc + phscat identical (discrete exact, energy to rtol)
        @test !isempty(emit_full)
        @test length(emit_full) == length(emit_sing)
        @test all(((p, q),) -> ceq(p, q), zip(emit_full, emit_sing))   # same LORs through the core
    end

    @testset "coincidence HDF5 round-trip" begin
        dir = mktempdir(); p = joinpath(dir, "lors.h5")
        w = CoincidenceWriter(p, Dict("has_randoms" => false); block=32)   # small block → many flushes
        ref = Tuple[]
        for ev in 1:150
            args = (ev, 100.0 + ev, -50.0, 30.0, 511.0, 0.0, ev % 20, ev % 48,
                    -200.0, 40.0, -30.0, 480.0, 0.0, (ev + 3) % 20, (ev + 5) % 48,
                    0.25 * ev, 1.0, 2.0, 3.0, Int8(ev % 2))   # dt_ns at index 16, then x0,y0,z0,truth
            push_coincidence!(w, args...); push!(ref, args)
        end
        @test close(w) == 150
        got = Tuple[]
        foreach_coincidences_hdf5(p; batch=32) do b
            for i in 1:length(b)
                push!(got, (Int(b.event[i]), decode_xyz(b.x1[i]), decode_e(b.e1[i]),
                            Int(b.iz1[i]), Int(b.iphi1[i]), Int(b.iz2[i]), Int(b.iphi2[i]),
                            Int8(b.truth[i]), Float64(b.dt[i])))
            end
        end
        @test length(got) == 150
        ok = true
        for (r, g) in zip(ref, got)
            ok &= g[1] == r[1] && g[8] == r[20]                        # event, truth exact
            ok &= g[4] == r[7] && g[5] == r[8] && g[6] == r[14] && g[7] == r[15]  # blocks exact
            ok &= abs(g[2] - r[2]) <= XYZ_SCALE_MM/2 + 1e-6             # x1 within ½ step
            ok &= abs(g[3] - r[5]) <= E_SCALE_KEV/2 + 1e-6             # e1 within ½ step
            ok &= abs(g[9] - r[16]) <= 1e-3                            # dt_ns round-trips (Float32)
        end
        @test ok
        rm(dir; recursive=true)
    end

    @testset "activity model — truncated-exponential event times" begin
        m = ActivityModel(; t0=0.0, t1=600.0, half_life_s=O15_HALFLIFE_S, seed=7)
        T = m.t1 - m.t0; λ = m.λ
        N = 200_000
        ts = [event_time(m, ev) for ev in 1:N]

        @test all(t -> m.t0 <= t <= m.t1, ts)                      # within the window
        @test event_time(m, 12345) == event_time(m, 12345)        # deterministic
        @test event_time(ActivityModel(; seed=7), 99) == event_time(m, 99)   # reproducible per (seed, ev)
        @test event_time(ActivityModel(; seed=8), 99) != event_time(m, 99)   # seed changes it

        # Sampled mean ≈ the analytic truncated-exponential mean (over [t0,t1], rate λ).
        mean_analytic = m.t0 + 1/λ - T*exp(-λ*T) / (1 - exp(-λ*T))
        @test isapprox(sum(ts)/N, mean_analytic; atol=2.0)        # ~3 s of slack over 2e5 draws

        # Front-loaded: the fraction in the first half matches the analytic CDF at the midpoint.
        cdf_mid = (1 - exp(-λ*T/2)) / (1 - exp(-λ*T))
        frac_first = count(t -> t <= m.t0 + T/2, ts) / N
        @test frac_first > 0.5 && isapprox(frac_first, cdf_mid; atol=0.01)

        event_time(m, 1)                                          # compile
        @test (@allocated event_time(m, 7)) == 0                 # stateless → no heap

        @test ActivityModel(Dict("timing" => Dict("t1_s" => 300.0))).t1 == 300.0
        @test ActivityModel(Dict{String,Any}()).t1 == 600.0      # default 10 min, ¹⁵O
        @test_throws ErrorException ActivityModel(; t0=10.0, t1=5.0)
    end

    @testset "scintillation timing — first-photon jitter + TOF + stamp" begin
        csi = load_material(DATA_DIR, "CsI")            # single component, 1000 ns
        bgo = load_material(DATA_DIR, "BGO")            # two components, 747 / 2253 ns
        w   = load_material(DATA_DIR, "Water")          # non-scintillator
        rng = MersenneTwister(2024)
        E = 0.511                                       # MeV

        # Mean jitter ≈ 1 / (N_det · r0), with r0 = Σ wₖ/τₖ over the decay mixture.
        N = 200_000
        for mat in (csi, bgo)
            r0 = sum(mat.scint_decay_w ./ mat.scint_decay_ns)
            N_det = mat.light_yield * E * mat.pde
            mean_analytic = 1.0 / (N_det * r0)
            js = [first_photon_jitter(mat, E, rng) for _ in 1:N]
            @test all(>=(0.0), js)
            @test isapprox(sum(js)/N, mean_analytic; rtol=0.02)
        end
        # BGO r0 is the explicit weighted sum (cross-check the mixture handling).
        @test isapprox(sum(bgo.scint_decay_w ./ bgo.scint_decay_ns),
                       0.34/747.0 + 0.66/2253.0; rtol=1e-12)

        # Non-scintillator (no yield / PDE) → exactly zero jitter.
        @test first_photon_jitter(w, E, rng) == 0.0

        # TOF: a known separation over c.
        @test isapprox(tof_ns((0.0,0.0,0.0), (C_MM_PER_NS,0.0,0.0)), 1.0; rtol=1e-12)
        @test tof_ns((1.0,2.0,3.0), (1.0,2.0,3.0)) == 0.0

        # Full stamp = t_annih (s→ns) + TOF + jitter; with zero jitter (Water) it's exact.
        ts = photon_timestamp(2.0, (0.0,0.0,0.0), (C_MM_PER_NS,0.0,0.0), E, w, rng)
        @test isapprox(ts, 2.0e9 + 1.0; rtol=1e-15)

        # Allocation-free hot path.
        first_photon_jitter(csi, E, rng); photon_timestamp(0.0, (0.,0.,0.), (1.,0.,0.), E, csi, rng)
        @test (@allocated first_photon_jitter(csi, E, rng)) == 0
        @test (@allocated photon_timestamp(0.0, (0.,0.,0.), (1.,0.,0.), E, csi, rng)) == 0
    end

    @testset "randoms pairing — cross-event window" begin
        # Singles on the absolute clock [ns]; pair cross-event ones within τ=1.5, skip same-event.
        t_abs = [0.0, 1.0, 2.0, 100.0, 100.5, 100.8]
        ev    = Int32[1,   1,   2,   3,     4,     3]      # (1,2) same decay; (4,6) same decay
        pairs = Tuple{Int,Int,Float64}[]
        n = pair_randoms((i, j, Δ) -> push!(pairs, (i, j, Δ)), t_abs, ev, 1.5)

        # Expected: (2,3) Δ1.0  [t=1↔2];  (4,5) Δ0.5  [100↔100.5];  (5,6) Δ0.3  [100.5↔100.8].
        # NOT (1,2): same decay.  NOT (1,3): Δ=2 > τ.  NOT (4,6): same decay (skipped, window continues).
        @test n == 3
        @test (2, 3, 1.0) in pairs
        @test any(p -> p[1]==4 && p[2]==5 && isapprox(p[3], 0.5), pairs)
        @test any(p -> p[1]==5 && p[2]==6 && isapprox(p[3], 0.3), pairs)
        @test all(((i,j,_),) -> ev[i] != ev[j], pairs)    # never same-event
        @test all(((_,_,Δ),) -> 0 < Δ <= 1.5, pairs)      # within the window

        # τ=0 → no pairs; a wide τ pairs every cross-event pair (here all but the 2 same-event ones).
        @test pair_randoms((i,j,Δ)->nothing, t_abs, ev, 0.0) == 0
        # Unsorted input is handled (internal sortperm): same result on a shuffled copy.
        perm = [4, 1, 5, 2, 6, 3]
        m = pair_randoms((i,j,Δ)->nothing, t_abs[perm], ev[perm], 1.5)
        @test m == 3
    end

    @testset "uniform volume sampling" begin
        rng = MersenneTwister(7)

        # Every drawn point lies inside its solid, for each shape.
        cyl = Cylinder(8.0, 8.0); sph = Sphere(8.0); box = Box(2.4, 2.4, 1.85)
        all_in = true
        for _ in 1:2000
            all_in &= is_inside(cyl, sample_point_in(cyl, rng))
            all_in &= is_inside(sph, sample_point_in(sph, rng))
            all_in &= is_inside(box, sample_point_in(box, rng))
        end
        @test all_in

        # Uniform-by-volume signatures: a uniform ball has mean radius 3R/4; a uniform
        # cylinder is centred on its axis with <z> ≈ 0.
        N = 40000
        rsum = 0.0; zsum = 0.0; xsum = 0.0; ysum = 0.0
        for _ in 1:N
            p = sample_point_in(sph, rng)
            rsum += sqrt(p[1]^2 + p[2]^2 + p[3]^2)
            q = sample_point_in(cyl, rng)
            xsum += q[1]; ysum += q[2]; zsum += q[3]
        end
        @test isapprox(rsum / N, 0.75 * 8.0; rtol=0.02)     # <r> = 3R/4 = 6.0
        @test abs(zsum / N) < 0.15 && abs(xsum / N) < 0.15 && abs(ysum / N) < 0.15

        # Through a placed PhysicalVolume the source samples in the world frame.
        pv = PhysicalVolume(LogicalVolume("ph", sph, load_material(DATA_DIR, "Water")),
                            (10.0, 0.0, 0.0))
        src = UniformVolumeSource(pv)
        @test all(is_inside(pv, sample_position(src, rng)) for _ in 1:2000)
    end

    @testset "back-to-back emission + acollinearity" begin
        rng = MersenneTwister(9)
        pv  = PhysicalVolume(LogicalVolume("ph", Sphere(8.0), load_material(DATA_DIR, "Water")),
                             (0.0, 0.0, 0.0))
        src = UniformVolumeSource(pv)

        # A point source returns its fixed position, with the same emission machinery.
        ps = PointSource((1.0, 2.0, 3.0))
        p0, _, _ = emit_pair(ps, rng; acol_fwhm_deg=0.0)
        @test p0 == (1.0, 2.0, 3.0)

        # Regression: the acollinearity near-pole path. rotate_to_global with a tuple
        # local vector and an axis along ±z must return a unit Vector{Float64}, not a tuple
        # (a tuple tripped the ::Vector{Float64} return when a photon was emitted ≈ along z).
        for axis in ((0.0, 0.0, 1.0), (0.0, 0.0, -1.0))
            g = PTCryspMC.rotate_to_global((0.02, -0.01, sqrt(1 - 0.02^2 - 0.01^2)), axis)
            @test g isa Vector{Float64}
            @test isapprox(sum(abs2, g), 1.0; atol=1e-9)
        end

        # Exactly back to back when acollinearity is off (accumulate, assert once).
        inside = true; unit = true; antiparallel = true
        for _ in 1:500
            pos, d1, d2 = emit_pair(src, rng; acol_fwhm_deg=0.0)
            inside &= is_inside(pv, pos)
            unit &= isapprox(d1[1]^2 + d1[2]^2 + d1[3]^2, 1.0) &&
                    isapprox(d2[1]^2 + d2[2]^2 + d2[3]^2, 1.0)
            antiparallel &= isapprox(d1[1]*d2[1] + d1[2]*d2[2] + d1[3]*d2[3], -1.0; atol=1e-9)
        end
        @test inside
        @test unit
        @test antiparallel

        # With 0.5° FWHM the deviation of d2 from −d1 has the right RMS angle; in the same
        # loop, the emission points fill the sphere phantom uniformly (the driver's source
        # path: emit_pair → sample_position), so their mean radius is 3R/4.
        fwhm = 0.5
        N = 20000; sumsq = 0.0; rsum = 0.0
        for _ in 1:N
            pos, d1, d2 = emit_pair(src, rng; acol_fwhm_deg=fwhm)
            rsum += sqrt(pos[1]^2 + pos[2]^2 + pos[3]^2)
            dot = -(d1[1]*d2[1] + d1[2]*d2[2] + d1[3]*d2[3])         # ≈ cos(deviation)
            δ = acos(clamp(dot, -1.0, 1.0))                         # deviation from 180°
            sumsq += δ^2
        end
        rms = sqrt(sumsq / N)
        σ = deg2rad(fwhm) / 2.3548200450309493
        @test isapprox(rms, sqrt(2) * σ; rtol=0.05)        # 2-D Gaussian: <δ²> = 2σ²
        @test isapprox(rsum / N, 0.75 * 8.0; rtol=0.02)    # uniform fill: <r> = 3R/4 = 6.0

        # emit_pair is heap-free per call — the acollinearity tilt routes through the tuple twin
        # rotate_to_global_t (not the Vector-returning rotate_to_global), so the source draw doesn't
        # reintroduce GC into the production MT loop. Per-annihilation allocation = the slope = 0.
        _emit_alloc_run(src, MersenneTwister(9), 10)        # compile
        b1000 = @allocated _emit_alloc_run(src, MersenneTwister(9), 1000)
        b2000 = @allocated _emit_alloc_run(src, MersenneTwister(9), 2000)
        @test (b2000 - b1000) == 0
    end

    @testset "detector response (smearing)" begin
        rng = MersenneTwister(13)

        # Energy resolution: FWHM(E) = a·√(511/E) → equals a (fractional) at 511 keV and
        # scales as √(511/E). a = 5% → 25.55 keV FWHM at 511, ×√2 wider at 255.5.
        @test isapprox(energy_fwhm(511.0, 0.05), 0.05 * 511.0)              # 25.55 keV
        @test isapprox(energy_fwhm(255.5, 0.05) / 255.5, 0.05 * sqrt(2.0)) # fractional ×√2
        @test isapprox(energy_sigma(511.0, 0.05), 0.05 * 511.0 / 2.3548200450309493)

        # smear_energy: a = 0 is exact; a > 0 is unbiased with the right σ.
        @test smear_energy(511.0, 0.0, rng) == 511.0
        N = 60000; s = 0.0; s2 = 0.0
        for _ in 1:N
            e = smear_energy(511.0, 0.05, rng); s += e; s2 += (e - 511.0)^2
        end
        @test isapprox(s / N, 511.0; atol=0.5)                              # unbiased
        @test isapprox(sqrt(s2 / N), energy_sigma(511.0, 0.05); rtol=0.03)  # right width

        # smear_position: σ = 0 is exact; σ > 0 is unbiased per axis with width σ.
        @test smear_position((1.0, 2.0, 3.0), 0.0, rng) === (1.0, 2.0, 3.0)
        sx = 0.0; sx2 = 0.0
        for _ in 1:N
            p = smear_position((10.0, 0.0, -5.0), 1.7, rng); sx += p[1]; sx2 += (p[1] - 10.0)^2
        end
        @test isapprox(sx / N, 10.0; atol=0.05)
        @test isapprox(sqrt(sx2 / N), 1.7; rtol=0.03)
    end

    @testset "run config (TOML)" begin
        path = tempname() * ".toml"
        write(path, """
        [source]
        nevents = 1234
        [output]
        tag = "myrun"
        """)
        cfg = read_config(path)
        @test cfg_get(cfg, "source", "nevents", 0) == 1234
        @test cfg_get(cfg, "source", "missing", 99) == 99    # default for a missing key
        @test cfg_get(cfg, "nope", "x", 7) == 7              # default for a missing section
        @test run_tag(cfg, path) == "myrun"                  # [output].tag wins
        rm(path)

        # Without [output].tag, the tag is the config filename (without .toml).
        p2 = joinpath(tempdir(), "sphere_water_csi.toml")
        write(p2, "[source]\nnevents = 1\n")
        @test run_tag(read_config(p2), p2) == "sphere_water_csi"
        rm(p2)
    end
end
