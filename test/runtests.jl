using PTCryspMC
using Test
using Random

const DATA_DIR = joinpath(@__DIR__, "..", "data")
const GEOM_JSON = joinpath(@__DIR__, "..", "geometry", "geometry.json")

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
        bad = tmpgeo("""{"shape":"sphere","radius_cm":5.0,"material":"Water"}""")
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
end
