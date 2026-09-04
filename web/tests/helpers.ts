import type {
  AtomRef,
  EntityRole,
  PreparedEntityFeatures,
  PreparedSite,
  ResidueRef,
  Vec3,
} from "../src/index.js";

let nextIndex = 1;

export function resetAtomIndices(): void {
  nextIndex = 1;
}

export function residue(entity: EntityRole, overrides: Partial<ResidueRef> = {}): ResidueRef {
  return {
    name: entity === "protein" ? "ALA" : entity === "ligand" ? "LIG" : "HOH",
    number: entity === "protein" ? 10 : entity === "ligand" ? 1 : 100,
    chain: entity === "protein" ? "A" : entity === "ligand" ? "L" : "W",
    ...overrides,
  };
}

export function atom(
  entity: EntityRole,
  position: Vec3,
  overrides: Partial<AtomRef> = {},
): AtomRef {
  const index = nextIndex++;
  return {
    index,
    parentIndex: index,
    name: entity === "water" ? "O" : "C",
    atomicNumber: entity === "water" ? 8 : 6,
    position,
    residue: residue(entity),
    entity,
    ...overrides,
  };
}

export function emptyFeatures(): PreparedEntityFeatures {
  return {
    hydrophobic: [],
    acceptors: [],
    donors: [],
    rings: [],
    charges: [],
    halogenAcceptors: [],
    halogenDonors: [],
    metalTargets: [],
  };
}

export function emptySite(): PreparedSite {
  return {
    protein: emptyFeatures(),
    ligand: emptyFeatures(),
    metals: [],
    waters: [],
  };
}
