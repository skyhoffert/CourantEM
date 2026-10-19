# 0008 — Rule-based mesher and KiCad import

Status: implemented (`web/src/mesh/mesher.ts`, `web/src/import/kicad.ts`, `web/src/import/boardsim.ts`).

## Mesher

Per axis, from features:

- **Sheet edges** (zero-thickness copper) with the side the metal is on → thirds-rule line pair (note 0002).
- **Planes** (dielectric interfaces, board edges, thick-copper faces) → a line exactly there.
- **Spans** that need ≥ n cells (CPW gaps, dielectric layers).
- **Max-cell regions** (λ/20 at fmax in each dielectric, via barrels).
- **Refinement points** where the cell is `hmin`.
- **Manual extra lines.**

These become a size field: every constraint is a cap that relaxes away from where it applies with slope 0.8·(ratio − 1), so the field is graded by construction. Each interval between fixed lines gets `n = ceil(∫dx/s)` cells placed at equal steps of that integral, so no cell exceeds the field and no splitting pass is needed. (A first version split oversized cells afterwards; it cascaded and blew up.) On the GCPW cross-section: max neighbour ratio 1.25 (y), 1.27 (z), all rules met. `estimate()` reports cells, Δt, steps and memory before a run.

## KiCad import

`importKicad` reads the `.kicad_pcb` S-expression (KiCad 6–8): stackup (copper and core/prepreg layers with thickness, εr, tan δ), nets, track segments, vias, filled zone polygons, footprint pads, and the Edge.Cuts bounding box. Constructs from real boards (`web/test/kicad.test.ts`): arc tracks (split at ≤10°), pads as rect, round-rect (corner radius from `roundrect_rratio`), circle, oval (a stadium), custom (anchor plus `gr_poly` primitives, as their convex hull), trapezoid/chamfered as rect; plated through-hole pads also get a barrel like a via; the outline takes lines, rectangles, arcs, circles and polygons on Edge.Cuts. Coordinates become metres with y up; z is stacked bottom-up from the bottom copper.

`buildBoardSim` turns the board into a simulation:

- mesh from copper edges (axis-aligned track sides and zone edges, with metal side from point-in-polygon), board edges, dielectric planes, via positions;
- dielectrics per layer (Djordjevic–Sarkar fit when lossy), copper loss when lossy;
- copper rasterised per layer onto Yee edges (an edge is metal if both ends and the midpoint are on copper; tracks are exact stadium shapes), row by row: each mesh row's copper is computed once as x-intervals (even-odd polygon crossings, stadium chords) that reproduce the point test exactly, so the cost is rows × objects rather than samples × objects;
- via barrels as solid cells inside the drill radius between their layers;
- the bottom copper must be a solid plane (it becomes the PEC floor) — milestone-1 scope;
- mode ports where a net's track leaves the board's left/right edge. The board is extruded past those edges by a feed length so source and monitor sit on a uniform line; via rows are carried through the feed as walls, otherwise the feed's top grounds would float. Only vias of the bottom plane's net within 2 mm of the port plane become walls, so stitching vias elsewhere or signal vias never short the feed.

## Gerber + drill import, stackup editor

`importGerber` (`web/src/import/gerber.ts`) builds the same `Board` model from RS-274X copper layers, a profile layer and Excellon drill files, so everything downstream (mesher, rasteriser, ports, 3D view) is shared with the KiCad path.

- **Gerber:** `%FS` (leading/trailing zero omission), `%MO`, standard apertures C/R/O/P, aperture macros (primitives 1, 4, 5, 20, 21 with `$n` arithmetic — enough for KiCad's RoundRect), D01/D02/D03, linear and circular interpolation (G74/G75, arcs linearised at 7.5°), regions G36/G37. Draws with a circular aperture become tracks (stadiums, exactly as KiCad segments); flashes and regions become polygons; draws with other apertures become the convex hull of the swept shape. Clear polarity, step-and-repeat and aperture transforms are reported as warnings, not applied.
- **Excellon:** METRIC/INCH with LZ/TZ and explicit `000.000` formats, decimal coordinates, tool table, `; TYPE=PLATED|NON_PLATED`. Slots (G85) are skipped.
- **File roles** come from X2 `%TF.FileFunction` when present, otherwise from KiCad (`-F_Cu`, `-In1_Cu`, `-Edge_Cuts`), Altium (`.GTL/.GBL/.Gn/.GKO/.GM1`) and generic names. Copper is ordered top to bottom.
- **Vias:** every plated hole is a via barrel; pad flashes centred on a hole become its lands (largest diameter = land size) and are not rasterised as separate copper. A hole without lands is a through via.
- **Nets:** X2 `%TO.N` names when present. Otherwise connectivity: objects on one layer are joined when they touch (sweep over x with bounding boxes, then exact stadium/polygon tests), vias join every layer between their lands, and components are named `Net-n`.
- **Outline:** bounding box of the profile layer; the copper bounding box (with a warning) if there is none.

Gerbers carry no stackup, so the app has a **stackup editor**: a table, top to bottom, of copper rows (thickness, rms roughness Rq) and dielectric rows (thickness, εr, tan δ at 10 GHz). It starts from the KiCad file's stackup or, for Gerbers, an evenly split 1.6 mm FR-4 board, and every edit re-imports the board. Roughness feeds the Groisse model of the copper SIBC (one model per board, the roughest copper wins — the SIBC is global for now).

A top solder mask can be added in the editor (KiCad's `Top Solder Mask` stackup layer is read automatically). It is modelled as a uniform dielectric slab from the top copper surface up by its thickness, filling the gaps; with thick copper the bars override it, so the mask on top of the copper is not modelled. On the test GCPW a 20 µm, εr 3.5 mask lowers Z0 by ~1.1 Ω; the same slab in the Laplace reference agrees within 0.8% (`web/test/thick.test.ts`).

Test (`web/test/gerber.test.ts`): parser units (formats, units, apertures, the RoundRect macro, an arc region's area, Excellon LZ/TZ) and a round trip — the KiCad GCPW fixture written out as KiCad-style Gerbers, with and without X2 net attributes, imports back to the same copper on a 48k-point raster, the same mesh edges and vias, one ground net, and a simulation whose first 300 steps are bit-identical to the KiCad import's. The browser test uploads the microstrip fixture as Gerbers, enters its stackup in the editor and checks the solved line.

## Simulated band

For large boards, a band of ± W around the port tracks can be simulated instead of the whole board (copper and vias outside dropped). Where the band cuts through the board, the cross-section continues into the y CPML instead of ending in air: a first version ended the band's pours and substrate in an air margin, and the parallel-plate region between the top pours and the bottom plane then formed an open-edged, nearly lossless cavity — the realistic coupon below rang at −20 dB for nanoseconds. Extending the cut edges into the absorber lets those waves leave (−58 dB by 410 ps), and on the GCPW fixture a ±1.5 mm band gives the same port Z0 as the whole board to 0.001 Ω.

`web/test/launch.test.ts` runs a realistic coupon (`test/fixtures/launch.kicad_pcb`: edge launches with round-rect and plated through-hole pins in rotated footprints, a jog of 45° segments and an arc, stepped-back top pour, via fences and stitching vias, solder mask, rounded outline) through the whole path: port planes 3.5 mm in, ±3 mm band, stackup thickness and mask. It decays to −51 dB in 361 ps, Z0 46.5 Ω, |S11| ≤ −21.5 dB, |S21| −0.25 dB at 20 GHz.

## Port planes

Each port has a plane, by default at the board edge, that the user can move into the board (port editor: "Left/Right port plane, mm in"). The board is cut there: copper, vias and mesh lines beyond the plane are dropped, the cross-section at the plane is extruded through the feed and into the CPML, and S-parameters are referred to the plane. This is how connector launches are excluded, matching TRL-calibrated measurements whose reference planes sit on the line. The port's net must have an x-directed track crossing the plane; that track sets the mode-port current loop.

**Full 2-port.** By default a run excites port 1 and the export assumes a symmetric, reciprocal line (S22 = S11, S12 = S21). The "full 2-port" option runs a second simulation exciting port 2 for boards that are not symmetric. On a one-sided pad (asymmetric, radiating open microstrip) the two runs give S12 = S21 to 1e-4 at 2 GHz and ~6e-3 at 40 GHz: space and surface waves from the radiating pad reach the monitor planes (the error drops with longer feeds and is not the far-end termination — solving S = B·A⁻¹ from both runs' four waves gives the same). The phases of S11 and S22 differ by up to 177°, as they should.

Test (`web/test/ports.test.ts`): the microstrip fixture with a wide pad 0.2–0.8 mm from each edge. Planes at 1 mm: |S11| ≤ −50 dB and S21 within 0.0000 dB / 0.045° of the plain line with the same planes. Planes at the edges: the pads show as |S11| up to −1.6 dB.

## Results

| Fixture | Result |
|---|---|
| Microstrip, 0.55 mm on 10-mil RO4350B, 8 mm | |S11| ≤ −49 dB, |S21| within 0.012 dB, Z0 49.90 Ω vs Hammerstad–Jensen 50.38 Ω (−1.0%) at 40 µm cells |
| GCPW with via fences (0.3 mm drills, 0.8 mm pitch), 10 mm | |S11| ≤ −26 dB, |S21| within 0.16 dB, Z0 49.7 Ω |

The GCPW-with-vias mismatch and loss are physical: the feeds carry ideal walls while the board has a discrete fence, and some energy leaks into the parallel-plate region outside the fence.

## Limits (next)

Ports only on the left/right sides (planes may sit inside the board) and on x-directed tracks; Gerber clear polarity and step-and-repeat not applied; solder mask only as a slab (not conformal over copper); arcs and angled tracks are staircased; board outline taken as its bounding box; no thick copper from the importer (stackup thickness is read but copper is modelled as sheets).
