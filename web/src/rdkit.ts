import type { JSMol, RDKitModule } from "@rdkit/rdkit";

import { angle, vector } from "./geometry.js";
import type {
  AromaticRing,
  AtomRef,
  ChargeGroup,
  EntityRole,
  HalogenAcceptor,
  HalogenDonor,
  HydrogenBondAcceptor,
  HydrogenBondDonor,
  HydrophobicAtom,
  MetalAtom,
  MetalTarget,
  PreparedEntityFeatures,
  PreparedSite,
  Vec3,
  WaterMolecule,
} from "./types.js";

export interface RdkitComponentInput {
  /** Sanitized RDKit.js molecule whose atom order matches `atoms`. */
  readonly molecule: JSMol;
  /** Original structure identities and 3D coordinates in RDKit atom order. */
  readonly atoms: readonly AtomRef[];
}

export interface RdkitSiteInput {
  readonly protein: readonly RdkitComponentInput[];
  readonly ligand: readonly RdkitComponentInput[];
  readonly waters?: readonly AtomRef[];
}

interface RdkitJsonAtom {
  readonly z?: number;
  readonly impHs?: number;
  readonly chg?: number;
}

interface RdkitJsonBond {
  readonly atoms?: readonly number[];
  readonly bo?: number;
}

interface RdkitRepresentation {
  readonly name?: string;
  readonly formatVersion?: number;
  readonly aromaticAtoms?: readonly number[];
  readonly atomRings?: readonly (readonly number[])[];
}

interface RdkitJson {
  readonly defaults?: {
    readonly atom?: RdkitJsonAtom;
    readonly bond?: RdkitJsonBond;
  };
  readonly molecules?: readonly {
    readonly atoms?: readonly RdkitJsonAtom[];
    readonly bonds?: readonly RdkitJsonBond[];
    readonly extensions?: readonly RdkitRepresentation[];
  }[];
}

interface ComponentTopology {
  readonly atoms: readonly AtomRef[];
  readonly atomicNumbers: readonly number[];
  readonly implicitHydrogens: readonly number[];
  readonly formalCharges: readonly number[];
  readonly neighbors: readonly (readonly number[])[];
  readonly aromaticAtoms: ReadonlySet<number>;
  readonly atomRings: readonly (readonly number[])[];
  readonly sp3Atoms: ReadonlySet<number>;
}

const ACCEPTOR_SMARTS = "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=[O,N,P,S])]),$([nH0,o,s;+0])]";
const DONOR_SMARTS = "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),$([n;H1;+0])]";
const METAL_ATOMIC_NUMBERS = new Set([
  3, 11, 12, 13, 19, 20, 24, 25, 26, 27, 28, 29, 30, 31, 37, 38, 44, 45, 46, 47,
  48, 49, 51, 55, 56, 57, 58, 59, 62, 63, 64, 65, 70, 71, 74, 76, 77, 78, 79,
  80, 81, 82,
]);
const ELEMENT_SYMBOLS = new Map<number, string>([
  [3, "LI"], [11, "NA"], [12, "MG"], [13, "AL"], [19, "K"], [20, "CA"],
  [24, "CR"], [25, "MN"], [26, "FE"], [27, "CO"], [28, "NI"], [29, "CU"],
  [30, "ZN"], [31, "GA"], [37, "RB"], [38, "SR"], [44, "RU"], [45, "RH"],
  [46, "PD"], [47, "AG"], [48, "CD"], [49, "IN"], [51, "SB"], [55, "CS"],
  [56, "BA"], [57, "LA"], [58, "CE"], [59, "PR"], [62, "SM"], [63, "EU"],
  [64, "GD"], [65, "TB"], [70, "YB"], [71, "LU"], [74, "W"], [76, "OS"],
  [77, "IR"], [78, "PT"], [79, "AU"], [80, "HG"], [81, "TL"], [82, "PB"],
]);
const HALOGEN_SYMBOLS = new Map<number, string>([
  [9, "F"], [17, "CL"], [35, "BR"], [53, "I"],
]);
const AROMATIC_AMINO_ACIDS = new Set(["TYR", "TRP", "HIS", "PHE"]);
const SIDECHAIN_METAL_ATOMS = new Map<string, ReadonlySet<string>>([
  ["ASP", new Set(["OD1", "OD2"])],
  ["GLU", new Set(["OE1", "OE2"])],
  ["SER", new Set(["OG"])],
  ["THR", new Set(["OG1"])],
  ["TYR", new Set(["OH"])],
  ["HIS", new Set(["ND1", "NE2"])],
  ["CYS", new Set(["SG"])],
]);
const POSITIVE_PROTEIN_ATOMS = new Map<string, ReadonlySet<string>>([
  ["ARG", new Set(["NE", "NH1", "NH2"])],
  ["HIS", new Set(["ND1", "NE2"])],
  ["LYS", new Set(["NZ"])],
]);
const NEGATIVE_PROTEIN_ATOMS = new Map<string, ReadonlySet<string>>([
  ["ASP", new Set(["OD1", "OD2"])],
  ["GLU", new Set(["OE1", "OE2"])],
]);

function parseJson(payload: string): RdkitJson {
  let parsed: unknown;
  try {
    parsed = JSON.parse(payload);
  } catch (error) {
    throw new Error("RDKit returned invalid molecular JSON", { cause: error });
  }
  if (typeof parsed !== "object" || parsed === null) {
    throw new TypeError("RDKit molecular JSON must be an object");
  }
  return parsed as RdkitJson;
}

function parseMatches(payload: string, atomCount: number, smarts: string): readonly (readonly number[])[] {
  let parsed: unknown;
  try {
    parsed = JSON.parse(payload);
  } catch (error) {
    throw new Error(`RDKit returned invalid match JSON for ${smarts}`, { cause: error });
  }
  if (typeof parsed === "object" && parsed !== null && !Array.isArray(parsed)
      && Object.keys(parsed).length === 0) return [];
  if (!Array.isArray(parsed)) throw new TypeError(`RDKit match result for ${smarts} is not an array`);
  return parsed.map((match, matchIndex) => {
    const indices = (match as { atoms?: unknown }).atoms;
    if (!Array.isArray(indices) || indices.length === 0) {
      throw new TypeError(`RDKit match ${matchIndex} for ${smarts} has no atoms`);
    }
    return indices.map((index) => {
      if (!Number.isInteger(index) || (index as number) < 0 || (index as number) >= atomCount) {
        throw new RangeError(`RDKit match ${matchIndex} for ${smarts} has invalid atom ${String(index)}`);
      }
      return index as number;
    });
  });
}

function matches(rdkit: RDKitModule, molecule: JSMol, atomCount: number, smarts: string): readonly (readonly number[])[] {
  const query = rdkit.get_qmol(smarts);
  if (query === null) throw new SyntaxError(`RDKit could not parse SMARTS: ${smarts}`);
  try {
    return parseMatches(molecule.get_substruct_matches(query), atomCount, smarts);
  } finally {
    query.delete();
  }
}

function topology(rdkit: RDKitModule, component: RdkitComponentInput): ComponentTopology {
  const json = parseJson(component.molecule.get_json());
  const molecule = json.molecules?.[0];
  const rawAtoms = molecule?.atoms;
  if (rawAtoms === undefined || rawAtoms.length !== component.atoms.length) {
    throw new RangeError(`RDKit molecule has ${rawAtoms?.length ?? 0} atoms but ${component.atoms.length} identities were supplied`);
  }
  const atomDefaults = json.defaults?.atom ?? {};
  const atomicNumbers = rawAtoms.map((atom, index) => atom.z ?? atomDefaults.z ?? component.atoms[index]!.atomicNumber);
  atomicNumbers.forEach((atomicNumber, index) => {
    if (atomicNumber !== component.atoms[index]!.atomicNumber) {
      throw new RangeError(
        `RDKit atom ${index} has atomic number ${atomicNumber} but the supplied identity has ${component.atoms[index]!.atomicNumber}`,
      );
    }
  });
  const implicitHydrogens = rawAtoms.map((atom) => atom.impHs ?? atomDefaults.impHs ?? 0);
  const formalCharges = rawAtoms.map((atom) => atom.chg ?? atomDefaults.chg ?? 0);
  const adjacency = Array.from({ length: rawAtoms.length }, () => [] as number[]);
  for (const [bondIndex, rawBond] of (molecule?.bonds ?? []).entries()) {
    const endpoints = rawBond.atoms;
    if (!Array.isArray(endpoints) || endpoints.length !== 2) {
      throw new TypeError(`RDKit bond ${bondIndex} does not have two endpoints`);
    }
    const [left, right] = endpoints;
    if (!Number.isInteger(left) || !Number.isInteger(right)
        || left! < 0 || right! < 0 || left! >= rawAtoms.length || right! >= rawAtoms.length) {
      throw new RangeError(`RDKit bond ${bondIndex} has invalid endpoints`);
    }
    adjacency[left!]!.push(right!);
    adjacency[right!]!.push(left!);
  }
  const representation = molecule?.extensions?.find((extension) => extension.name === "rdkitRepresentation")
    ?? molecule?.extensions?.find((extension) => extension.atomRings !== undefined);
  const aromaticAtoms = new Set([
    ...(representation?.aromaticAtoms ?? []),
    ...matches(rdkit, component.molecule, rawAtoms.length, "[a]").map((match) => match[0]!),
  ]);
  const atomRings = (representation?.atomRings ?? []).filter((ring) =>
    ring.every((index) => Number.isInteger(index) && index >= 0 && index < rawAtoms.length));
  const sp3Atoms = new Set(matches(rdkit, component.molecule, rawAtoms.length, "[*^3]").map((match) => match[0]!));
  return {
    atoms: component.atoms,
    atomicNumbers,
    implicitHydrogens,
    formalCharges,
    neighbors: adjacency,
    aromaticAtoms,
    atomRings,
    sp3Atoms,
  };
}

function centroid(atoms: readonly AtomRef[]): Vec3 {
  if (atoms.length === 0) throw new RangeError("Cannot calculate the centroid of an empty atom set");
  const sum = atoms.reduce<Vec3>((value, atom) => [
    value[0] + atom.position[0], value[1] + atom.position[1], value[2] + atom.position[2],
  ], [0, 0, 0]);
  return [sum[0] / atoms.length, sum[1] / atoms.length, sum[2] / atoms.length];
}

function cross(left: Vec3, right: Vec3): Vec3 {
  return [
    left[1] * right[2] - left[2] * right[1],
    left[2] * right[0] - left[0] * right[2],
    left[0] * right[1] - left[1] * right[0],
  ];
}

function normalized(value: Vec3): Vec3 | undefined {
  const length = Math.hypot(...value);
  if (length === 0) return undefined;
  return [value[0] / length, value[1] / length, value[2] / length];
}

function ringNormal(atoms: readonly AtomRef[]): Vec3 | undefined {
  if (atoms.length < 5) return undefined;
  return normalized(cross(
    vector(atoms[0]!.position, atoms[2]!.position),
    vector(atoms[4]!.position, atoms[0]!.position),
  ));
}

function ringIsPlanar(indices: readonly number[], data: ComponentTopology): boolean {
  const ringSet = new Set(indices);
  const normals: Vec3[] = [];
  for (const index of indices) {
    const ringNeighbors = data.neighbors[index]!.filter((neighbor) => ringSet.has(neighbor));
    if (ringNeighbors.length !== 2) return false;
    const first = vector(data.atoms[index]!.position, data.atoms[ringNeighbors[0]!]!.position);
    const second = vector(data.atoms[index]!.position, data.atoms[ringNeighbors[1]!]!.position);
    const normal = normalized(cross(first, second));
    if (normal === undefined) return false;
    normals.push(normal);
  }
  for (const left of normals) {
    for (const right of normals) {
      const separation = angle(left, right);
      if (separation > 5 && separation < 175) return false;
    }
  }
  return true;
}

function groupKey(atom: AtomRef): string {
  const insertionCode = atom.residue.insertionCode ?? "";
  return `${atom.residue.chain}\u0000${atom.residue.number}\u0000${insertionCode}\u0000${atom.residue.name}`;
}

function groupAtoms(atoms: readonly AtomRef[]): Map<string, AtomRef[]> {
  const groups = new Map<string, AtomRef[]>();
  for (const atom of atoms) {
    const key = groupKey(atom);
    const group = groups.get(key);
    if (group === undefined) groups.set(key, [atom]);
    else group.push(atom);
  }
  return groups;
}

function proteinCharges(atoms: readonly AtomRef[]): ChargeGroup[] {
  const charges: ChargeGroup[] = [];
  for (const residueAtoms of groupAtoms(atoms).values()) {
    const residueName = residueAtoms[0]!.residue.name.toUpperCase();
    for (const [charge, names] of [
      ["positive", POSITIVE_PROTEIN_ATOMS.get(residueName)],
      ["negative", NEGATIVE_PROTEIN_ATOMS.get(residueName)],
    ] as const) {
      if (names === undefined) continue;
      const contributors = residueAtoms.filter((atom) => names.has(atom.name.trim().toUpperCase()));
      if (contributors.length > 0) charges.push({
        atoms: contributors,
        center: centroid(contributors),
        charge,
        functionalGroup: residueName,
      });
    }
  }
  return charges;
}

function ligandCharges(data: ComponentTopology): ChargeGroup[] {
  const output: ChargeGroup[] = [];
  const seen = new Set<string>();
  const add = (indices: readonly number[], charge: ChargeGroup["charge"], functionalGroup: string, centerIndices = indices): void => {
    const sorted = [...new Set(indices)].sort((left, right) => left - right);
    const key = `${charge}:${functionalGroup}:${sorted.join(",")}`;
    if (seen.has(key) || sorted.length === 0) return;
    seen.add(key);
    const atoms = sorted.map((index) => data.atoms[index]!);
    const centerAtoms = centerIndices.map((index) => data.atoms[index]!);
    const group: ChargeGroup = { atoms, center: centroid(centerAtoms), charge, functionalGroup };
    if (functionalGroup === "tertamine") {
      const neighbors = data.neighbors[sorted[0]!]!.slice(0, 3).map((index) => data.atoms[index]!);
      const normal = neighbors.length === 3 ? normalized(cross(
        vector(neighbors[0]!.position, neighbors[1]!.position),
        vector(neighbors[2]!.position, neighbors[0]!.position),
      )) : undefined;
      if (normal !== undefined) output.push({ ...group, normal });
      else output.push(group);
    } else output.push(group);
  };

  data.atoms.forEach((_atom, index) => {
    const atomicNumber = data.atomicNumbers[index]!;
    const neighbors = data.neighbors[index]!;
    const neighborElements = neighbors.map((neighbor) => data.atomicNumbers[neighbor]!);
    if (atomicNumber === 7 && neighborElements.every((number) => number !== 1) && neighbors.length === 4) {
      add([index], "positive", "quartamine");
    } else if (atomicNumber === 7 && data.sp3Atoms.has(index) && neighbors.length + data.implicitHydrogens[index]! >= 3) {
      add([index], "positive", "tertamine");
    }
    if (atomicNumber === 16 && neighborElements.every((number) => number !== 1) && neighbors.length === 3) {
      add([index], "positive", "sulfonium");
    }
    if (atomicNumber === 15 && neighbors.length > 0 && neighborElements.every((number) => number === 8)) {
      add([index, ...neighbors], "negative", "phosphate", [index]);
    }
    if (atomicNumber === 16 && neighborElements.filter((number) => number === 8).length === 3) {
      add([index, ...neighbors.filter((neighbor) => data.atomicNumbers[neighbor] === 8)], "negative", "sulfonicacid", [index]);
    } else if (atomicNumber === 16 && neighborElements.filter((number) => number === 8).length === 4) {
      add([index, ...neighbors], "negative", "sulfate", [index]);
    }
    if (atomicNumber === 6 && neighborElements.filter((number) => number === 8).length === 2
        && neighborElements.filter((number) => number === 6).length === 1) {
      const oxygens = neighbors.filter((neighbor) => data.atomicNumbers[neighbor] === 8);
      add(oxygens, "negative", "carboxylate");
    } else if (atomicNumber === 6 && neighbors.length === 3 && neighborElements.every((number) => number === 7)
        && Math.min(...neighbors.map((neighbor) => data.neighbors[neighbor]!.length)) === 1) {
      add(neighbors, "positive", "guanidine", [index]);
    }
  });
  return output;
}

function metalTargets(data: ComponentTopology, role: "protein" | "ligand"): MetalTarget[] {
  const output: MetalTarget[] = [];
  const seen = new Set<number>();
  const add = (index: number, type: string, location: MetalTarget["location"]): void => {
    if (seen.has(index)) return;
    seen.add(index);
    output.push({ atom: data.atoms[index]!, type, location });
  };
  data.atoms.forEach((atom, index) => {
    const element = data.atomicNumbers[index]!;
    if (role === "protein") {
      const residue = atom.residue.name.toUpperCase();
      const name = atom.name.trim().toUpperCase();
      if (SIDECHAIN_METAL_ATOMS.get(residue)?.has(name)) {
        add(index, ELEMENT_SYMBOLS.get(element) ?? (element === 8 ? "O" : element === 7 ? "N" : "S"), "protein");
      } else if (name === "O" && residue !== "HOH") add(index, "O", "protein");
      return;
    }
    const neighbors = data.neighbors[index]!;
    const neighborElements = neighbors.map((neighbor) => data.atomicNumbers[neighbor]!);
    if (element === 8) {
      if ((neighborElements.filter((number) => number === 1).length === 1 || data.implicitHydrogens[index] === 1)
          && neighbors.length + data.implicitHydrogens[index]! === 2) add(index, "O", "ligand");
      if (neighbors.some((neighbor) => data.aromaticAtoms.has(neighbor)) && !data.aromaticAtoms.has(index)) {
        add(index, "O", "ligand");
      }
    }
    if (element === 6 && neighborElements.filter((number) => number === 8).length === 2
        && neighborElements.filter((number) => number === 6).length === 1) {
      for (const neighbor of neighbors.filter((neighbor) => data.atomicNumbers[neighbor] === 8)) add(neighbor, "O", "ligand");
    }
    if (element === 15 && (neighborElements.filter((number) => number === 8).length === 2
        || neighborElements.filter((number) => number === 8).length >= 3)) {
      for (const neighbor of neighbors.filter((neighbor) => data.atomicNumbers[neighbor] === 8)) add(neighbor, "O", "ligand");
    }
    if (element === 7 && neighborElements.filter((number) => number === 6).length === 2) add(index, "N", "ligand");
    if (element === 16 && neighbors.some((neighbor) => data.aromaticAtoms.has(neighbor)) && !data.aromaticAtoms.has(index)) {
      add(index, "S", "ligand");
    }
    if (element === 16 && neighborElements.length > 0 && neighborElements.every((number) => number === 26)) {
      add(index, "S", "ligand");
    }
  });
  return output;
}

function prepareComponent(rdkit: RDKitModule, component: RdkitComponentInput, role: "protein" | "ligand") {
  const data = topology(rdkit, component);
  const acceptorIndices = new Set(matches(rdkit, component.molecule, data.atoms.length, ACCEPTOR_SMARTS).map((match) => match[0]!));
  const donorIndices = new Set(matches(rdkit, component.molecule, data.atoms.length, DONOR_SMARTS).map((match) => match[0]!));
  const hydrophobic: HydrophobicAtom[] = [];
  const acceptors: HydrogenBondAcceptor[] = [];
  const donors: HydrogenBondDonor[] = [];
  const rings: AromaticRing[] = [];
  const halogenAcceptors: HalogenAcceptor[] = [];
  const halogenDonors: HalogenDonor[] = [];

  data.atoms.forEach((atom, index) => {
    const element = data.atomicNumbers[index]!;
    const neighbors = data.neighbors[index]!;
    const neighborElements = neighbors.map((neighbor) => data.atomicNumbers[neighbor]!);
    if (element === 6 && neighborElements.every((number) => number === 1 || number === 6)) {
      hydrophobic.push({ atom, neighbors: neighbors.map((neighbor) => data.atoms[neighbor]!.index) });
    }
    if (acceptorIndices.has(index) && ![9, 17, 35, 53].includes(element)) acceptors.push({ atom, type: "regular" });
    if (donorIndices.has(index)) {
      for (const neighbor of neighbors.filter((neighbor) => data.atomicNumbers[neighbor] === 1)) {
        donors.push({ atom, hydrogen: data.atoms[neighbor]!, type: "regular" });
      }
    }
    if ([7, 8, 16].includes(element)) {
      const proximal = neighbors.filter((neighbor) => [6, 7, 15, 16].includes(data.atomicNumbers[neighbor]!));
      if (proximal.length === 1) halogenAcceptors.push({
        oxygen: atom,
        neighbor: data.atoms[proximal[0]!]!,
        type: "regular",
      });
    }
    if ([9, 17, 35, 53].includes(element)) {
      const carbons = neighbors.filter((neighbor) => data.atomicNumbers[neighbor] === 6);
      if (carbons.length === 1) halogenDonors.push({
        halogen: atom,
        carbon: data.atoms[carbons[0]!]!,
        type: HALOGEN_SYMBOLS.get(element)!,
      });
    }
  });

  for (const [ringIndex, indices] of data.atomRings.entries()) {
    if (indices.length <= 4 || indices.length > 6) continue;
    const atoms = indices.map((index) => data.atoms[index]!);
    const aromatic = indices.every((index) => data.aromaticAtoms.has(index));
    const aminoAcid = AROMATIC_AMINO_ACIDS.has(atoms[0]!.residue.name.toUpperCase());
    if (!aromatic && !aminoAcid && !ringIsPlanar(indices, data)) continue;
    const normal = ringNormal(atoms);
    if (normal === undefined) continue;
    rings.push({ id: `${groupKey(atoms[0]!)}:${ringIndex}`, atoms, center: centroid(atoms), normal });
  }

  return {
    data,
    features: {
      hydrophobic,
      acceptors,
      donors,
      rings,
      charges: role === "protein" ? proteinCharges(data.atoms) : ligandCharges(data),
      halogenAcceptors,
      halogenDonors,
      metalTargets: metalTargets(data, role),
    } satisfies PreparedEntityFeatures,
  };
}

/**
 * Prepare PLIP-compatible features from sanitized RDKit.js components.
 *
 * RDKit supplies topology, aromaticity, valence, and donor/acceptor perception;
 * the original structure remains authoritative for atom identity and 3D
 * coordinates. Hydrogen-bond events require explicit hydrogen atoms in both
 * the RDKit molecule and the supplied atom array, matching PLIP's geometric
 * donor-angle contract.
 */
export function prepareRdkitSite(rdkit: RDKitModule, input: RdkitSiteInput): PreparedSite {
  const proteinComponents = input.protein.map((component) => prepareComponent(rdkit, component, "protein"));
  const ligandComponents = input.ligand.map((component) => prepareComponent(rdkit, component, "ligand"));
  const merge = (components: readonly ReturnType<typeof prepareComponent>[]): PreparedEntityFeatures => ({
    hydrophobic: components.flatMap(({ features }) => features.hydrophobic),
    acceptors: components.flatMap(({ features }) => features.acceptors),
    donors: components.flatMap(({ features }) => features.donors),
    rings: components.flatMap(({ features }) => features.rings),
    charges: components.flatMap(({ features }) => features.charges),
    halogenAcceptors: components.flatMap(({ features }) => features.halogenAcceptors),
    halogenDonors: components.flatMap(({ features }) => features.halogenDonors),
    metalTargets: components.flatMap(({ features }) => features.metalTargets),
  });
  const waterAtoms = (input.waters ?? []).filter((atom) => atom.atomicNumber === 8);
  const waters: WaterMolecule[] = waterAtoms.map((oxygen) => ({ oxygen }));
  const protein = merge(proteinComponents);
  const ligand = merge(ligandComponents);
  const waterTargets: MetalTarget[] = waterAtoms.map((atom) => ({ atom, type: "O", location: "water" }));
  const metals: MetalAtom[] = [...proteinComponents, ...ligandComponents].flatMap(({ data }) =>
    data.atoms.flatMap((atom, index) => METAL_ATOMIC_NUMBERS.has(data.atomicNumbers[index]!) ? [{
      atom,
      element: ELEMENT_SYMBOLS.get(data.atomicNumbers[index]!) ?? String(data.atomicNumbers[index]),
    }] : []));
  return {
    protein,
    ligand: { ...ligand, metalTargets: [...ligand.metalTargets, ...waterTargets] },
    metals,
    waters,
  };
}

/** Useful when a caller needs the same feature preparation for one side only. */
export function prepareRdkitEntity(
  rdkit: RDKitModule,
  components: readonly RdkitComponentInput[],
  role: Exclude<EntityRole, "water">,
): PreparedEntityFeatures {
  const prepared = components.map((component) => prepareComponent(rdkit, component, role));
  return {
    hydrophobic: prepared.flatMap(({ features }) => features.hydrophobic),
    acceptors: prepared.flatMap(({ features }) => features.acceptors),
    donors: prepared.flatMap(({ features }) => features.donors),
    rings: prepared.flatMap(({ features }) => features.rings),
    charges: prepared.flatMap(({ features }) => features.charges),
    halogenAcceptors: prepared.flatMap(({ features }) => features.halogenAcceptors),
    halogenDonors: prepared.flatMap(({ features }) => features.halogenDonors),
    metalTargets: prepared.flatMap(({ features }) => features.metalTargets),
  };
}
