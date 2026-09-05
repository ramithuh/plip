import { beforeEach, describe, expect, it } from "vitest";
import {
  detectCationPi,
  detectHalogenBonds,
  detectHydrogenBonds,
  detectHydrophobic,
  detectMetalComplexes,
  detectPiStacking,
  detectPreparedSite,
  detectSaltBridges,
  detectWaterBridges,
  PLIP_DEFAULTS,
  PlipFingerprint,
  type AromaticRing,
  type ChargeGroup,
  type HydrogenBondAcceptor,
  type HydrogenBondDonor,
} from "../src/index.js";
import { atom, emptyFeatures, emptySite, resetAtomIndices, residue } from "./helpers.js";

beforeEach(resetAtomIndices);

describe("prepared-feature interaction detectors", () => {
  it("uses PLIP's strict hydrophobic distance boundaries and patch refinement", () => {
    const protein = atom("protein", [0, 0, 0]);
    const ligandNear = atom("ligand", [3, 0, 0]);
    const ligandFar = atom("ligand", [3.5, 0, 0]);
    const ligandDisconnected = atom("ligand", [0, 3.2, 0]);
    const events = detectHydrophobic(
      [{ atom: protein, neighbors: [] }],
      [
        { atom: ligandNear, neighbors: [ligandFar.index] },
        { atom: ligandFar, neighbors: [ligandNear.index] },
        { atom: ligandDisconnected, neighbors: [] },
      ],
    );
    expect(events.map((event) => event.ligandAtoms[0]!.index).sort()).toEqual([
      ligandNear.index,
      // Native PLIP's edge-only clustering omits the disconnected contact.
    ].sort());
    expect(detectHydrophobic(
      [{ atom: protein, neighbors: [] }],
      [{ atom: atom("ligand", [4, 0, 0]), neighbors: [] }],
    )).toHaveLength(0);
  });

  it("keeps the largest-angle hydrogen bond for each donor", () => {
    const donor = atom("ligand", [0, 0, 0], { name: "N", atomicNumber: 7 });
    const hydrogen = atom("ligand", [1, 0, 0], { name: "H", atomicNumber: 1 });
    const straight = atom("protein", [2.8, 0, 0], { name: "O", atomicNumber: 8 });
    const bent = atom("protein", [2.8, 0.6, 0], { name: "O", atomicNumber: 8 });
    const donors: HydrogenBondDonor[] = [{ atom: donor, hydrogen, type: "N" }];
    const acceptors: HydrogenBondAcceptor[] = [
      { atom: bent, type: "O" },
      { atom: straight, type: "O" },
    ];
    const events = detectHydrogenBonds(acceptors, [], [], donors);
    expect(events).toHaveLength(1);
    expect(events[0]!.direction).toBe("ligand-donor");
    expect(events[0]!.proteinAtoms[0]!.index).toBe(straight.index);
    expect(events[0]!.geometry?.donorAngle).toBeCloseTo(180);
  });

  it("detects parallel and perpendicular pi stacking and rejects exact limits", () => {
    const proteinAtom = atom("protein", [0, 0, 0], { residue: residue("protein", { name: "PHE" }) });
    const ligandAtom = atom("ligand", [0, 0, 4]);
    const parallel: AromaticRing = {
      id: "l-parallel",
      atoms: [ligandAtom],
      center: [0, 0, 4],
      normal: [Math.sin(Math.PI / 18), 0, Math.cos(Math.PI / 18)],
    };
    const perpendicular: AromaticRing = {
      id: "l-perpendicular",
      atoms: [ligandAtom],
      center: [0, 0, 4],
      normal: [1, 0, 0],
    };
    const proteinRing: AromaticRing = {
      id: "p",
      atoms: [proteinAtom],
      center: [0, 0, 0],
      normal: [0, 0, 1],
    };
    expect(detectPiStacking([proteinRing], [parallel])[0]!.subtype).toBe("parallel");
    expect(detectPiStacking([proteinRing], [perpendicular])[0]!.subtype).toBe("perpendicular");
    expect(detectPiStacking([proteinRing], [{ ...parallel, center: [0, 0, 5.5] }])).toHaveLength(0);
  });

  it("detects cation-pi contacts and enforces tertiary-amine orientation", () => {
    const ringAtom = atom("protein", [0, 0, 0], { residue: residue("protein", { name: "TRP" }) });
    const cationAtom = atom("ligand", [0, 0, 4], { name: "N", atomicNumber: 7 });
    const ring: AromaticRing = { id: "p", atoms: [ringAtom], center: [0, 0, 0], normal: [0, 0, 1] };
    const charge: ChargeGroup = {
      atoms: [cationAtom], center: [0, 0, 4], charge: "positive", functionalGroup: "tertamine", normal: [0, 0, 1],
    };
    const events = detectCationPi([ring], [], [], [charge]);
    expect(events).toHaveLength(1);
    expect(events[0]!.direction).toBe("ligand-cation");
    expect(detectCationPi([ring], [], [], [{ ...charge, normal: [1, 0, 0] }])).toHaveLength(0);
  });

  it("detects only opposite-charge salt bridges", () => {
    const proteinAtom = atom("protein", [0, 0, 0], { name: "NZ", atomicNumber: 7 });
    const ligandAtom = atom("ligand", [3, 0, 0], { name: "O", atomicNumber: 8 });
    const positive: ChargeGroup = { atoms: [proteinAtom], center: proteinAtom.position, charge: "positive" };
    const negative: ChargeGroup = { atoms: [ligandAtom], center: ligandAtom.position, charge: "negative" };
    expect(detectSaltBridges([positive], [negative])).toHaveLength(1);
    expect(detectSaltBridges([positive], [{ ...negative, charge: "positive" }])).toHaveLength(0);
  });

  it("enforces PLIP halogen distance and angular geometry", () => {
    const oxygen = atom("protein", [0, 0, 0], { name: "O", atomicNumber: 8 });
    const neighbor = atom("protein", [-0.5, Math.sqrt(3) / 2, 0]);
    const halogen = atom("ligand", [3, 0, 0], { name: "CL", atomicNumber: 17 });
    const carbon = atom("ligand", [3 + Math.cos(Math.PI / 12), Math.sin(Math.PI / 12), 0]);
    expect(detectHalogenBonds(
      [{ oxygen, neighbor, type: "O" }],
      [{ halogen, carbon, type: "halocarbon" }],
    )).toHaveLength(1);
    expect(detectHalogenBonds(
      [{ oxygen, neighbor, type: "O" }],
      [{ halogen: { ...halogen, position: [4, 0, 0] }, carbon, type: "halocarbon" }],
    )).toHaveLength(0);
  });

  it("retains at most two water bridges and removes direct-donor duplicates", () => {
    const waterO = atom("water", [0, 0, 0]);
    const makeAcceptor = (number: number, angleDegrees: number): HydrogenBondAcceptor => {
      const radians = angleDegrees * Math.PI / 180;
      return {
        atom: atom("ligand", [3 * Math.cos(radians), 3 * Math.sin(radians), 0], {
          name: "O", atomicNumber: 8, residue: residue("ligand", { number }),
        }),
        type: "O",
      };
    };
    const hydrogen = atom("protein", [1, 0, 0], { name: "H", atomicNumber: 1 });
    const donorAtom = atom("protein", [3, 0, 0], { name: "N", atomicNumber: 7 });
    const donor: HydrogenBondDonor = { atom: donorAtom, hydrogen, type: "N" };
    const acceptors = [makeAcceptor(1, 110), makeAcceptor(2, 105), makeAcceptor(3, 130)];
    const events = detectWaterBridges([], [donor], acceptors, [], [{ oxygen: waterO }]);
    expect(events).toHaveLength(2);
    expect(events.map((event) => event.ligandResidue.number).sort()).toEqual([1, 2]);
    expect(detectWaterBridges([], [donor], acceptors, [], [{ oxygen: waterO }], PLIP_DEFAULTS, [{
      family: "HydrogenBond",
      proteinResidue: donorAtom.residue,
      ligandResidue: acceptors[0]!.atom.residue,
      proteinAtoms: [donorAtom],
      ligandAtoms: [acceptors[0]!.atom],
      distance: 3,
      direction: "protein-donor",
    }])).toHaveLength(0);
  });

  it("reports protein-ligand metal contacts but not water-only or ligand-internal contacts", () => {
    const metal = atom("ligand", [0, 0, 0], { name: "FE", atomicNumber: 26 });
    const proteinTarget = atom("protein", [2.2, 0, 0], { name: "NE2", atomicNumber: 7 });
    const waterTarget = atom("water", [0, 2.2, 0]);
    const ligandTarget = atom("ligand", [0, 0, 2.2], { name: "O", atomicNumber: 8 });
    const events = detectMetalComplexes(
      [{ atom: metal, element: "FE" }],
      [
        { atom: proteinTarget, type: "N", location: "protein" },
        { atom: waterTarget, type: "O", location: "water" },
        { atom: ligandTarget, type: "O", location: "ligand" },
      ],
    );
    expect(events).toHaveLength(1);
    expect(events[0]!.proteinAtoms[0]!.index).toBe(proteinTarget.index);
    expect(events[0]!.geometry?.observedCoordination).toBe(3);
  });

  it("reproduces PLIP's greedy metal-geometry choice", () => {
    const metal = atom("ligand", [0, 0, 0], { name: "ZN", atomicNumber: 30 });
    const tetrahedralPoints = [
      [1, 1, 1],
      [1, -1, -1],
      [-1, 1, -1],
      [-1, -1, 1],
    ] as const;
    const tetrahedral = detectMetalComplexes(
      [{ atom: metal, element: "ZN" }],
      tetrahedralPoints.map((position, index) => ({
        atom: atom("protein", position, {
          name: `T${index}`,
          residue: residue("protein", { number: 20 + index }),
        }),
        type: "N",
        location: "protein" as const,
      })),
    );
    expect(tetrahedral).toHaveLength(4);
    expect(tetrahedral[0]!.geometry).toMatchObject({
      coordination: 4,
      observedCoordination: 4,
      shape: "tetrahedral",
    });
    expect(tetrahedral[0]!.geometry?.rms).toBeLessThan(0.1);

    const squarePlanar = detectMetalComplexes(
      [{ atom: metal, element: "ZN" }],
      [[2, 0, 0], [-2, 0, 0], [0, 2, 0], [0, -2, 0]].map((position, index) => ({
        atom: atom("protein", position as [number, number, number], {
          name: `S${index}`,
          residue: residue("protein", { number: 30 + index }),
        }),
        type: "O",
        location: "protein" as const,
      })),
    );
    expect(squarePlanar).toHaveLength(4);
    expect(squarePlanar[0]!.geometry).toMatchObject({
      coordination: 4,
      observedCoordination: 4,
      // PLIP 3.0.1 itself selects tetrahedral for this ideal square-planar
      // input because its decision heuristic compares adjacent candidates.
      shape: "tetrahedral",
    });
    expect(squarePlanar[0]!.geometry?.rms).toBeCloseTo(75.70171728567325, 6);
  });
});

describe("PLIP reporting refinements", () => {
  it("suppresses a hydrogen bond that duplicates a salt bridge", () => {
    const proteinDonor = atom("protein", [0, 0, 0], { name: "NZ", atomicNumber: 7 });
    const hydrogen = atom("protein", [1, 0, 0], { name: "H", atomicNumber: 1 });
    const ligandAcceptor = atom("ligand", [2.8, 0, 0], { name: "O", atomicNumber: 8 });
    const site = emptySite();
    const prepared = {
      ...site,
      protein: {
        ...emptyFeatures(),
        donors: [{ atom: proteinDonor, hydrogen, type: "N" }],
        charges: [{ atoms: [proteinDonor], center: proteinDonor.position, charge: "positive" as const }],
      },
      ligand: {
        ...emptyFeatures(),
        acceptors: [{ atom: ligandAcceptor, type: "O" }],
        charges: [{ atoms: [ligandAcceptor], center: ligandAcceptor.position, charge: "negative" as const }],
      },
    };
    const events = detectPreparedSite(prepared);
    expect(events.map((event) => event.family)).toEqual(["SaltBridge"]);
  });

  it("suppresses ring hydrophobics and histidine cation-pi when stacking is present", () => {
    const proteinAtom = atom("protein", [0, 0, 0], { residue: residue("protein", { name: "HIS" }) });
    const ligandAtom = atom("ligand", [0, 0, 3.5]);
    const proteinRing: AromaticRing = {
      id: "p", atoms: [proteinAtom], center: [0, 0, 0], normal: [0, 0, 1],
    };
    const ligandRing: AromaticRing = {
      id: "l", atoms: [ligandAtom], center: [0, 0, 3.5], normal: [0.1, 0, 0.995],
    };
    const prepared = {
      ...emptySite(),
      protein: {
        ...emptyFeatures(),
        rings: [proteinRing],
        hydrophobic: [{ atom: proteinAtom, neighbors: [] }],
        charges: [{ atoms: [proteinAtom], center: proteinRing.center, charge: "positive" as const }],
      },
      ligand: {
        ...emptyFeatures(),
        rings: [ligandRing],
        hydrophobic: [{ atom: ligandAtom, neighbors: [] }],
      },
    };
    expect(detectPreparedSite(prepared).map((event) => event.family)).toEqual(["PiStacking"]);
  });

  it("exposes the core through the fingerprint API", () => {
    const protein = atom("protein", [0, 0, 0]);
    const ligand = atom("ligand", [3, 0, 0]);
    const result = new PlipFingerprint().detect({
      ...emptySite(),
      protein: { ...emptyFeatures(), hydrophobic: [{ atom: protein, neighbors: [] }] },
      ligand: { ...emptyFeatures(), hydrophobic: [{ atom: ligand, neighbors: [] }] },
    });
    expect(result.events).toHaveLength(1);
    expect(result.events[0]!.family).toBe("Hydrophobic");
  });
});
