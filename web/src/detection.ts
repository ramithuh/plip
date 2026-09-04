import { PLIP_DEFAULTS, type PlipThresholds } from "./config.js";
import { angle, distance, projectOntoPlane, vector } from "./geometry.js";
import type {
  AromaticRing,
  AtomRef,
  ChargeGroup,
  HalogenAcceptor,
  HalogenDonor,
  HydrogenBondAcceptor,
  HydrogenBondDonor,
  HydrophobicAtom,
  InteractionEvent,
  MetalAtom,
  MetalTarget,
  PreparedSite,
  ResidueRef,
  WaterMolecule,
} from "./types.js";

function belongsTo(atom: AtomRef, entity: "protein" | "ligand"): boolean {
  return atom.entity === entity;
}

function residues(proteinAtom: AtomRef, ligandAtom: AtomRef): {
  proteinResidue: ResidueRef;
  ligandResidue: ResidueRef;
} {
  return {
    proteinResidue: proteinAtom.residue,
    ligandResidue: ligandAtom.residue,
  };
}

function strictRange(value: number, minimum: number, maximum: number): boolean {
  return minimum < value && value < maximum;
}

function groupBy<T, K>(values: readonly T[], keyOf: (value: T) => K): Map<K, T[]> {
  const groups = new Map<K, T[]>();
  for (const value of values) {
    const key = keyOf(value);
    const group = groups.get(key);
    if (group === undefined) groups.set(key, [value]);
    else group.push(value);
  }
  return groups;
}

function atomSetEquals(left: readonly AtomRef[], right: readonly AtomRef[]): boolean {
  if (left.length !== right.length) return false;
  const indices = new Set(left.map((atom) => atom.parentIndex));
  return right.every((atom) => indices.has(atom.parentIndex));
}

export function detectHydrophobic(
  protein: readonly HydrophobicAtom[],
  ligand: readonly HydrophobicAtom[],
  thresholds: PlipThresholds = PLIP_DEFAULTS,
): InteractionEvent[] {
  const raw: Array<InteractionEvent & { protein: HydrophobicAtom; ligand: HydrophobicAtom }> = [];
  for (const proteinFeature of protein) {
    for (const ligandFeature of ligand) {
      const separation = distance(proteinFeature.atom.position, ligandFeature.atom.position);
      if (!strictRange(separation, thresholds.minimumDistance, thresholds.hydrophobicDistanceMax)) continue;
      raw.push({
        family: "Hydrophobic",
        ...residues(proteinFeature.atom, ligandFeature.atom),
        proteinAtoms: [proteinFeature.atom],
        ligandAtoms: [ligandFeature.atom],
        distance: separation,
        protein: proteinFeature,
        ligand: ligandFeature,
      });
    }
  }

  // PLIP first keeps the closest contact for each ligand atom/protein residue pair.
  const byLigandAndResidue = new Map<string, typeof raw[number]>();
  for (const event of raw) {
    const key = `${event.ligand.atom.index}|${event.protein.atom.residue.chain}|${event.protein.atom.residue.number}`;
    const previous = byLigandAndResidue.get(key);
    if (previous === undefined || event.distance < previous.distance) byLigandAndResidue.set(key, event);
  }

  // For one protein atom contacting a connected ligand patch, PLIP reports its shortest edge.
  const candidates = [...byLigandAndResidue.values()];
  const output: typeof candidates = [];
  const byProtein = groupBy(candidates, (event) => event.protein.atom.index);
  for (const events of byProtein.values()) {
    if (events.length === 1) {
      output.push(events[0]!);
      continue;
    }
    const eventByLigand = new Map(events.map((event) => [event.ligand.atom.index, event]));
    const visited = new Set<number>();
    for (const event of events) {
      if (visited.has(event.ligand.atom.index)) continue;
      const stack = [event.ligand.atom.index];
      const component: typeof events = [];
      while (stack.length > 0) {
        const atomIndex = stack.pop()!;
        if (visited.has(atomIndex)) continue;
        visited.add(atomIndex);
        const member = eventByLigand.get(atomIndex);
        if (member === undefined) continue;
        component.push(member);
        for (const neighbor of member.ligand.neighbors) {
          if (eventByLigand.has(neighbor)) stack.push(neighbor);
        }
      }
      output.push(component.reduce((best, candidate) => candidate.distance < best.distance ? candidate : best));
    }
  }
  return output.map(({ protein: _protein, ligand: _ligand, ...event }) => event);
}

function detectHydrogenBondDirection(
  acceptors: readonly HydrogenBondAcceptor[],
  donors: readonly HydrogenBondDonor[],
  proteinIsDonor: boolean,
  thresholds: PlipThresholds,
): InteractionEvent[] {
  const raw: InteractionEvent[] = [];
  for (const acceptor of acceptors) {
    for (const donor of donors) {
      const donorAcceptorDistance = distance(acceptor.atom.position, donor.atom.position);
      if (!strictRange(donorAcceptorDistance, thresholds.minimumDistance, thresholds.hydrogenBondDistanceMax)) continue;
      const donorAngle = angle(
        vector(donor.hydrogen.position, donor.atom.position),
        vector(donor.hydrogen.position, acceptor.atom.position),
      );
      if (!(donorAngle > thresholds.hydrogenBondDonorAngleMin)) continue;
      const proteinAtom = proteinIsDonor ? donor.atom : acceptor.atom;
      const ligandAtom = proteinIsDonor ? acceptor.atom : donor.atom;
      raw.push({
        family: "HydrogenBond",
        ...residues(proteinAtom, ligandAtom),
        proteinAtoms: [proteinAtom],
        ligandAtoms: [ligandAtom],
        distance: donorAcceptorDistance,
        direction: proteinIsDonor ? "protein-donor" : "ligand-donor",
        geometry: {
          distanceAH: distance(acceptor.atom.position, donor.hydrogen.position),
          donorAngle,
        },
      });
    }
  }

  // PLIP allows one hydrogen bond per donor and selects the largest donor angle.
  const byDonor = new Map<number, InteractionEvent>();
  for (const event of raw) {
    const donorIndex = (proteinIsDonor ? event.proteinAtoms[0] : event.ligandAtoms[0])!.index;
    const previous = byDonor.get(donorIndex);
    const currentAngle = event.geometry?.donorAngle as number;
    const previousAngle = previous?.geometry?.donorAngle as number | undefined;
    if (previous === undefined || currentAngle > previousAngle!) byDonor.set(donorIndex, event);
  }
  return [...byDonor.values()];
}

export function detectHydrogenBonds(
  proteinAcceptors: readonly HydrogenBondAcceptor[],
  proteinDonors: readonly HydrogenBondDonor[],
  ligandAcceptors: readonly HydrogenBondAcceptor[],
  ligandDonors: readonly HydrogenBondDonor[],
  thresholds: PlipThresholds = PLIP_DEFAULTS,
): InteractionEvent[] {
  return [
    ...detectHydrogenBondDirection(ligandAcceptors, proteinDonors, true, thresholds),
    ...detectHydrogenBondDirection(proteinAcceptors, ligandDonors, false, thresholds),
  ];
}

function refineHydrogenBondsAgainstSaltBridges(
  hydrogenBonds: readonly InteractionEvent[],
  saltBridges: readonly InteractionEvent[],
): InteractionEvent[] {
  return hydrogenBonds.filter((hydrogenBond) => {
    const proteinIndex = hydrogenBond.proteinAtoms[0]?.parentIndex;
    const ligandIndex = hydrogenBond.ligandAtoms[0]?.parentIndex;
    return !saltBridges.some((saltBridge) =>
      saltBridge.proteinAtoms.some((atom) => atom.parentIndex === proteinIndex)
      && saltBridge.ligandAtoms.some((atom) => atom.parentIndex === ligandIndex));
  });
}

export function detectPiStacking(
  proteinRings: readonly AromaticRing[],
  ligandRings: readonly AromaticRing[],
  thresholds: PlipThresholds = PLIP_DEFAULTS,
): InteractionEvent[] {
  const events: InteractionEvent[] = [];
  for (const proteinRing of proteinRings) {
    for (const ligandRing of ligandRings) {
      const separation = distance(proteinRing.center, ligandRing.center);
      if (!strictRange(separation, thresholds.minimumDistance, thresholds.piStackingDistanceMax)) continue;
      const rawAngle = angle(proteinRing.normal, ligandRing.normal);
      const ringAngle = Math.min(rawAngle, 180 - rawAngle);
      const proteinProjection = projectOntoPlane(ligandRing.normal, ligandRing.center, proteinRing.center);
      const ligandProjection = projectOntoPlane(proteinRing.normal, proteinRing.center, ligandRing.center);
      const offset = Math.min(
        distance(proteinProjection, ligandRing.center),
        distance(ligandProjection, proteinRing.center),
      );
      const parallel = 0 < ringAngle && ringAngle < thresholds.piStackingAngleDeviation;
      const perpendicular = 90 - thresholds.piStackingAngleDeviation < ringAngle
        && ringAngle < 90 + thresholds.piStackingAngleDeviation;
      if ((!parallel && !perpendicular) || !(offset < thresholds.piStackingOffsetMax)) continue;
      const proteinAtom = proteinRing.atoms[0]!;
      const ligandAtom = ligandRing.atoms[0]!;
      events.push({
        family: "PiStacking",
        ...residues(proteinAtom, ligandAtom),
        proteinAtoms: proteinRing.atoms,
        ligandAtoms: ligandRing.atoms,
        distance: separation,
        subtype: parallel ? "parallel" : "perpendicular",
        geometry: { angle: ringAngle, offset },
      });
    }
  }
  return events;
}

function detectCationPiDirection(
  rings: readonly AromaticRing[],
  charges: readonly ChargeGroup[],
  proteinIsCation: boolean,
  thresholds: PlipThresholds,
): InteractionEvent[] {
  const events: InteractionEvent[] = [];
  for (const ring of rings) {
    for (const charge of charges.filter((candidate) => candidate.charge === "positive")) {
      const separation = distance(ring.center, charge.center);
      const projection = projectOntoPlane(ring.normal, ring.center, charge.center);
      const offset = distance(projection, ring.center);
      if (!strictRange(separation, thresholds.minimumDistance, thresholds.cationPiDistanceMax)) continue;
      if (!(offset < thresholds.piStackingOffsetMax)) continue;
      if (charge.functionalGroup === "tertamine" && charge.normal !== undefined) {
        const rawAngle = angle(ring.normal, charge.normal);
        if (!(Math.min(rawAngle, 180 - rawAngle) < 30)) continue;
      }
      const ringAtom = ring.atoms[0]!;
      const chargeAtom = charge.atoms[0]!;
      const proteinAtom = proteinIsCation ? chargeAtom : ringAtom;
      const ligandAtom = proteinIsCation ? ringAtom : chargeAtom;
      events.push({
        family: "CationPi",
        ...residues(proteinAtom, ligandAtom),
        proteinAtoms: proteinIsCation ? charge.atoms : ring.atoms,
        ligandAtoms: proteinIsCation ? ring.atoms : charge.atoms,
        distance: separation,
        direction: proteinIsCation ? "protein-cation" : "ligand-cation",
        geometry: { offset },
      });
    }
  }
  return events;
}

export function detectCationPi(
  proteinRings: readonly AromaticRing[],
  proteinCharges: readonly ChargeGroup[],
  ligandRings: readonly AromaticRing[],
  ligandCharges: readonly ChargeGroup[],
  thresholds: PlipThresholds = PLIP_DEFAULTS,
): InteractionEvent[] {
  return [
    ...detectCationPiDirection(ligandRings, proteinCharges, true, thresholds),
    ...detectCationPiDirection(proteinRings, ligandCharges, false, thresholds),
  ];
}

function refineCationPiAgainstStacking(
  cationPi: readonly InteractionEvent[],
  stacking: readonly InteractionEvent[],
): InteractionEvent[] {
  return cationPi.filter((event) => {
    if (event.proteinResidue.name === "HIS" && event.direction === "protein-cation") {
      return !stacking.some((stack) =>
        stack.proteinResidue.chain === event.proteinResidue.chain
        && stack.proteinResidue.number === event.proteinResidue.number
        && atomSetEquals(stack.ligandAtoms, event.ligandAtoms));
    }
    if (event.ligandResidue.name === "HIS" && event.direction === "ligand-cation") {
      return !stacking.some((stack) =>
        stack.ligandResidue.chain === event.ligandResidue.chain
        && stack.ligandResidue.number === event.ligandResidue.number
        && atomSetEquals(stack.proteinAtoms, event.proteinAtoms));
    }
    return true;
  });
}

export function detectSaltBridges(
  proteinCharges: readonly ChargeGroup[],
  ligandCharges: readonly ChargeGroup[],
  thresholds: PlipThresholds = PLIP_DEFAULTS,
): InteractionEvent[] {
  const events: InteractionEvent[] = [];
  for (const proteinCharge of proteinCharges) {
    for (const ligandCharge of ligandCharges) {
      if (proteinCharge.charge === ligandCharge.charge) continue;
      const separation = distance(proteinCharge.center, ligandCharge.center);
      if (!strictRange(separation, thresholds.minimumDistance, thresholds.saltBridgeDistanceMax)) continue;
      const proteinAtom = proteinCharge.atoms[0]!;
      const ligandAtom = ligandCharge.atoms[0]!;
      events.push({
        family: "SaltBridge",
        ...residues(proteinAtom, ligandAtom),
        proteinAtoms: proteinCharge.atoms,
        ligandAtoms: ligandCharge.atoms,
        distance: separation,
        direction: proteinCharge.charge === "positive" ? "protein-cation" : "ligand-cation",
      });
    }
  }
  return events;
}

export function detectHalogenBonds(
  proteinAcceptors: readonly HalogenAcceptor[],
  ligandDonors: readonly HalogenDonor[],
  thresholds: PlipThresholds = PLIP_DEFAULTS,
): InteractionEvent[] {
  const events: InteractionEvent[] = [];
  for (const acceptor of proteinAcceptors) {
    for (const donor of ligandDonors) {
      const separation = distance(acceptor.oxygen.position, donor.halogen.position);
      if (!strictRange(separation, thresholds.minimumDistance, thresholds.halogenDistanceMax)) continue;
      const acceptorAngle = angle(
        vector(acceptor.oxygen.position, acceptor.neighbor.position),
        vector(acceptor.oxygen.position, donor.halogen.position),
      );
      const donorAngle = angle(
        vector(donor.halogen.position, acceptor.oxygen.position),
        vector(donor.halogen.position, donor.carbon.position),
      );
      if (!(thresholds.halogenAcceptorAngle - thresholds.halogenAngleDeviation < acceptorAngle
        && acceptorAngle < thresholds.halogenAcceptorAngle + thresholds.halogenAngleDeviation)) continue;
      if (!(thresholds.halogenDonorAngle - thresholds.halogenAngleDeviation < donorAngle
        && donorAngle < thresholds.halogenDonorAngle + thresholds.halogenAngleDeviation)) continue;
      events.push({
        family: "HalogenBond",
        ...residues(acceptor.oxygen, donor.halogen),
        proteinAtoms: [acceptor.oxygen],
        ligandAtoms: [donor.halogen],
        distance: separation,
        geometry: { acceptorAngle, donorAngle },
      });
    }
  }
  return events;
}

interface WaterBridgeCandidate {
  readonly acceptor: HydrogenBondAcceptor;
  readonly donor: HydrogenBondDonor;
  readonly water: WaterMolecule;
  readonly proteinIsDonor: boolean;
  readonly distanceAW: number;
  readonly distanceDW: number;
  readonly donorAngle: number;
  readonly waterAngle: number;
}

export function detectWaterBridges(
  proteinAcceptors: readonly HydrogenBondAcceptor[],
  proteinDonors: readonly HydrogenBondDonor[],
  ligandAcceptors: readonly HydrogenBondAcceptor[],
  ligandDonors: readonly HydrogenBondDonor[],
  waters: readonly WaterMolecule[],
  thresholds: PlipThresholds = PLIP_DEFAULTS,
  directHydrogenBonds: readonly InteractionEvent[] = [],
): InteractionEvent[] {
  const candidates: WaterBridgeCandidate[] = [];
  const directions: Array<{
    acceptors: readonly HydrogenBondAcceptor[];
    donors: readonly HydrogenBondDonor[];
    proteinIsDonor: boolean;
  }> = [
    { acceptors: ligandAcceptors, donors: proteinDonors, proteinIsDonor: true },
    { acceptors: proteinAcceptors, donors: ligandDonors, proteinIsDonor: false },
  ];
  for (const water of waters) {
    for (const direction of directions) {
      for (const acceptor of direction.acceptors) {
        const distanceAW = distance(acceptor.atom.position, water.oxygen.position);
        if (!(thresholds.waterBridgeDistanceMin <= distanceAW
          && distanceAW <= thresholds.waterBridgeDistanceMax)) continue;
        for (const donor of direction.donors) {
          const distanceDW = distance(donor.atom.position, water.oxygen.position);
          const donorAngle = angle(
            vector(donor.hydrogen.position, donor.atom.position),
            vector(donor.hydrogen.position, water.oxygen.position),
          );
          if (!(thresholds.waterBridgeDistanceMin <= distanceDW
            && distanceDW <= thresholds.waterBridgeDistanceMax
            && donorAngle > thresholds.waterBridgeThetaMin)) continue;
          const waterAngle = angle(
            vector(water.oxygen.position, acceptor.atom.position),
            vector(water.oxygen.position, donor.hydrogen.position),
          );
          if (!(thresholds.waterBridgeOmegaMin < waterAngle
            && waterAngle < thresholds.waterBridgeOmegaMax)) continue;
          candidates.push({
            acceptor,
            donor,
            water,
            proteinIsDonor: direction.proteinIsDonor,
            distanceAW,
            distanceDW,
            donorAngle,
            waterAngle,
          });
        }
      }
    }
  }

  // Retain one donor per water/acceptor and at most two bridges per water,
  // preferring water angles closest to PLIP's 110° optimum.
  const directDonorIndices = new Set(directHydrogenBonds.map((event) =>
    (event.direction === "protein-donor" ? event.proteinAtoms[0] : event.ligandAtoms[0])?.parentIndex));
  const byWaterAndAcceptor = new Map<string, WaterBridgeCandidate>();
  for (const candidate of candidates.filter((candidate) =>
    !directDonorIndices.has(candidate.donor.atom.parentIndex))) {
    const key = `${candidate.water.oxygen.index}|${candidate.acceptor.atom.index}`;
    const previous = byWaterAndAcceptor.get(key);
    if (previous === undefined
      || Math.abs(110 - candidate.waterAngle) < Math.abs(110 - previous.waterAngle)) {
      byWaterAndAcceptor.set(key, candidate);
    }
  }
  const byWater = groupBy([...byWaterAndAcceptor.values()], (candidate) => candidate.water.oxygen.index);
  return [...byWater.values()].flatMap((group) => [...group]
    .sort((left, right) => Math.abs(110 - left.waterAngle) - Math.abs(110 - right.waterAngle))
    .slice(0, 2)
    .map((candidate): InteractionEvent => {
      const proteinAtom = candidate.proteinIsDonor ? candidate.donor.atom : candidate.acceptor.atom;
      const ligandAtom = candidate.proteinIsDonor ? candidate.acceptor.atom : candidate.donor.atom;
      return {
        family: "WaterBridge",
        ...residues(proteinAtom, ligandAtom),
        proteinAtoms: [proteinAtom],
        ligandAtoms: [ligandAtom],
        distance: candidate.proteinIsDonor ? candidate.distanceAW : candidate.distanceDW,
        direction: candidate.proteinIsDonor ? "protein-donor" : "ligand-donor",
        geometry: {
          distanceAW: candidate.distanceAW,
          distanceDW: candidate.distanceDW,
          donorAngle: candidate.donorAngle,
          waterAngle: candidate.waterAngle,
        },
        water: candidate.water.oxygen,
      };
    }));
}

export function detectMetalComplexes(
  metals: readonly MetalAtom[],
  targets: readonly MetalTarget[],
  thresholds: PlipThresholds = PLIP_DEFAULTS,
): InteractionEvent[] {
  const events: InteractionEvent[] = [];
  for (const metal of metals) {
    const contacts = targets
      .map((target) => ({ target, separation: distance(metal.atom.position, target.atom.position) }))
      .filter(({ separation }) => separation < thresholds.metalDistanceMax);
    if (contacts.length === 0 || contacts.every(({ target }) => target.location === "water")) continue;
    for (const { target, separation } of contacts) {
      if (target.location === "water") continue;
      const targetIsProtein = belongsTo(target.atom, "protein");
      const metalIsProtein = belongsTo(metal.atom, "protein");
      if (!targetIsProtein && !metalIsProtein) continue;
      const proteinAtom = targetIsProtein ? target.atom : metal.atom;
      const ligandAtom = targetIsProtein ? metal.atom : target.atom;
      events.push({
        family: "MetalComplex",
        ...residues(proteinAtom, ligandAtom),
        proteinAtoms: [proteinAtom],
        ligandAtoms: [ligandAtom],
        distance: separation,
        subtype: target.type,
        geometry: { observedCoordination: contacts.length },
        metal: metal.atom,
      });
    }
  }
  return events;
}

export function detectPreparedSite(
  site: PreparedSite,
  thresholds: PlipThresholds = PLIP_DEFAULTS,
): InteractionEvent[] {
  const stacking = detectPiStacking(site.protein.rings, site.ligand.rings, thresholds);
  const saltBridges = detectSaltBridges(site.protein.charges, site.ligand.charges, thresholds);
  const hydrogenBonds = refineHydrogenBondsAgainstSaltBridges(detectHydrogenBonds(
    site.protein.acceptors,
    site.protein.donors,
    site.ligand.acceptors,
    site.ligand.donors,
    thresholds,
  ), saltBridges);
  const hydrophobic = detectHydrophobic(site.protein.hydrophobic, site.ligand.hydrophobic, thresholds)
    .filter((event) => !stacking.some((stack) =>
      stack.proteinAtoms.some((atom) => atom.index === event.proteinAtoms[0]?.index)
      && stack.ligandAtoms.some((atom) => atom.index === event.ligandAtoms[0]?.index)));
  const cationPi = refineCationPiAgainstStacking(detectCationPi(
    site.protein.rings,
    site.protein.charges,
    site.ligand.rings,
    site.ligand.charges,
    thresholds,
  ), stacking);
  return [
    ...hydrophobic,
    ...hydrogenBonds,
    ...stacking,
    ...cationPi,
    ...saltBridges,
    ...detectHalogenBonds(site.protein.halogenAcceptors, site.ligand.halogenDonors, thresholds),
    ...detectWaterBridges(
      site.protein.acceptors,
      site.protein.donors,
      site.ligand.acceptors,
      site.ligand.donors,
      site.waters,
      thresholds,
      hydrogenBonds,
    ),
    ...detectMetalComplexes(site.metals, [...site.ligand.metalTargets, ...site.protein.metalTargets], thresholds),
  ];
}
