export type Vec3 = readonly [x: number, y: number, z: number];
export type EntityRole = "protein" | "ligand" | "water";

export interface ResidueRef {
  readonly name: string;
  readonly number: number;
  readonly chain: string;
  readonly insertionCode?: string;
}

/** Atom identity and coordinates after structure preparation. */
export interface AtomRef {
  readonly index: number;
  readonly parentIndex: number;
  readonly name: string;
  readonly atomicNumber: number;
  readonly position: Vec3;
  readonly residue: ResidueRef;
  readonly entity: EntityRole;
}

export interface HydrophobicAtom {
  readonly atom: AtomRef;
  /** Directly bonded atom indices, used by PLIP's hydrophobic-contact refinement. */
  readonly neighbors: readonly number[];
}

export interface HydrogenBondAcceptor {
  readonly atom: AtomRef;
  readonly type: string;
}

export interface HydrogenBondDonor {
  readonly atom: AtomRef;
  readonly hydrogen: AtomRef;
  readonly type: string;
}

export interface AromaticRing {
  readonly id: string;
  readonly atoms: readonly AtomRef[];
  readonly center: Vec3;
  readonly normal: Vec3;
}

export interface ChargeGroup {
  readonly atoms: readonly AtomRef[];
  readonly center: Vec3;
  readonly charge: "positive" | "negative";
  readonly functionalGroup?: string;
  /** Plane normal for tertiary amines, when available. */
  readonly normal?: Vec3;
}

export interface HalogenAcceptor {
  readonly oxygen: AtomRef;
  readonly neighbor: AtomRef;
  readonly type: string;
}

export interface HalogenDonor {
  readonly halogen: AtomRef;
  readonly carbon: AtomRef;
  readonly type: string;
}

export interface MetalAtom {
  readonly atom: AtomRef;
  readonly element: string;
}

export interface MetalTarget {
  readonly atom: AtomRef;
  readonly type: string;
  readonly location: EntityRole;
}

export interface WaterMolecule {
  readonly oxygen: AtomRef;
}

export interface PreparedEntityFeatures {
  readonly hydrophobic: readonly HydrophobicAtom[];
  readonly acceptors: readonly HydrogenBondAcceptor[];
  readonly donors: readonly HydrogenBondDonor[];
  readonly rings: readonly AromaticRing[];
  readonly charges: readonly ChargeGroup[];
  readonly halogenAcceptors: readonly HalogenAcceptor[];
  readonly halogenDonors: readonly HalogenDonor[];
  readonly metalTargets: readonly MetalTarget[];
}

export interface PreparedSite {
  readonly protein: PreparedEntityFeatures;
  readonly ligand: PreparedEntityFeatures;
  readonly metals: readonly MetalAtom[];
  readonly waters: readonly WaterMolecule[];
}

export type InteractionFamily =
  | "Hydrophobic"
  | "HydrogenBond"
  | "PiStacking"
  | "CationPi"
  | "SaltBridge"
  | "HalogenBond"
  | "WaterBridge"
  | "MetalComplex";

export interface InteractionEvent {
  readonly family: InteractionFamily;
  readonly proteinResidue: ResidueRef;
  readonly ligandResidue: ResidueRef;
  readonly proteinAtoms: readonly AtomRef[];
  readonly ligandAtoms: readonly AtomRef[];
  readonly distance: number;
  readonly direction?: "protein-donor" | "ligand-donor" | "protein-cation" | "ligand-cation";
  readonly subtype?: string;
  readonly geometry?: Readonly<Record<string, number | string>>;
  readonly water?: AtomRef;
  readonly metal?: AtomRef;
}

export interface DetectionResult {
  readonly events: readonly InteractionEvent[];
}
