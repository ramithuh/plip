import { prepareTopologyComponent, type ComponentTopology } from "./rdkit.js";
import type { AtomRef, PreparedSite, Vec3 } from "./types.js";

/** JSON contract of the pinned original-C++ Open Babel WASM bridge. */
export interface OpenBabelAtom {
  readonly type?: string;
  readonly residueProperty8?: boolean;
  readonly index: number;
  readonly element: number;
  readonly formalCharge: number;
  readonly hybridization: number;
  readonly implicitHydrogens: number;
  readonly aromatic: boolean;
  readonly donor: boolean;
  readonly acceptor: boolean;
  readonly donorHydrogen: boolean;
  readonly xyz: Vec3;
}
export interface OpenBabelMolecule {
  readonly ok: true;
  readonly version: string;
  readonly atomsBefore: number;
  readonly atomsAfter: number;
  readonly addedHydrogens: number;
  readonly atoms: readonly OpenBabelAtom[];
  readonly bonds: readonly { readonly begin: number; readonly end: number; readonly order: number; readonly aromatic: boolean }[];
  readonly rings: readonly { readonly atoms: readonly number[]; readonly aromatic: boolean }[];
}

/**
 * Keep native OB perception; adapt only PLIP feature/identity representation.
 * undefined identities are contextual atoms, not selected receptor/ligand.
 * Shared functional-group/refinement rules still have fixture-bounded parity.
 */
export function prepareOpenBabelSite(
  molecule: OpenBabelMolecule,
  identities: readonly (AtomRef | undefined)[],
): PreparedSite {
  if (identities.length !== molecule.atoms.length) throw new RangeError("Open Babel identity count mismatch");
  const native = new Map(molecule.atoms.map((atom, i) => [atom.index, i]));
  if (native.size !== molecule.atoms.length) throw new RangeError("Duplicate Open Babel atom index");
  const adjacency = molecule.atoms.map(() => [] as number[]);
  for (const bond of molecule.bonds) {
    const a = native.get(bond.begin), b = native.get(bond.end);
    if (a === undefined || b === undefined) throw new RangeError("Invalid Open Babel bond endpoint");
    adjacency[a]!.push(b); adjacency[b]!.push(a);
  }
  identities.forEach((ref, i) => {
    if (ref && ref.atomicNumber !== molecule.atoms[i]!.element) throw new RangeError("Open Babel element/identity mismatch");
  });
  const prepare = (role: "protein" | "ligand") => {
    const selected = identities.flatMap((atom, i) => atom?.entity === role ? [i] : []);
    const local = new Map(selected.map((i, j) => [i, j]));
    const atoms = selected.map(i => identities[i]!);
    const raw = selected.map(i => molecule.atoms[i]!);
    const data: ComponentTopology = {
      atoms,
      atomicNumbers: raw.map(a => a.element),
      formalCharges: raw.map(a => a.formalCharge),
      implicitHydrogens: raw.map(a => a.implicitHydrogens),
      neighbors: selected.map(i => adjacency[i]!.flatMap(n => local.has(n) ? [local.get(n)!] : [])),
      aromaticAtoms: new Set(raw.flatMap((a, i) => a.aromatic ? [i] : [])),
      sp3Atoms: new Set(raw.flatMap((a, i) => a.hybridization === 3 ? [i] : [])),
      atomRings: molecule.rings.flatMap(r => {
        const indices = r.atoms.map(i => local.get(native.get(i)!));
        return indices.every(i => i !== undefined)
          ? [(indices as number[]).sort((a,b) => atoms[a]!.parentIndex - atoms[b]!.parentIndex)] : [];
      }),
    };
    const features = prepareTopologyComponent(data, role,
      new Set(raw.flatMap((a, i) => a.acceptor ? [i] : [])),
      new Set(raw.flatMap((a, i) => a.donor ? [i] : [])),
      new Set(raw.flatMap((a, i) => a.donorHydrogen ? [i] : []))).features;
    if (role !== 'protein' || raw.some(a => a.residueProperty8 === undefined)) return features;
    const groups = new Map<string, AtomRef[]>();
    atoms.forEach((atom,i) => {
      const name = atom.residue.name;
      const positive = ['ARG','HIS','LYS'].includes(name);
      const negative = ['ASP','GLU'].includes(name);
      if (!raw[i]!.residueProperty8 || !(positive && raw[i]!.type?.startsWith('N')
        || negative && raw[i]!.type?.startsWith('O'))) return;
      const key = `${atom.residue.chain}:${atom.residue.number}:${name}`;
      groups.set(key,[...(groups.get(key) ?? []),atom]);
    });
    const chargedNames = new Set(['ARG','HIS','LYS','ASP','GLU']);
    return {...features,charges:[...features.charges.filter(c => !chargedNames.has(c.atoms[0]!.residue.name)),
      ...[...groups.values()].map(atoms => ({atoms,center:atoms.reduce((s,a)=>
        s.map((v,i)=>v+a.position[i]!/atoms.length) as unknown as Vec3,[0,0,0] as Vec3),
        charge: (['ARG','HIS','LYS'].includes(atoms[0]!.residue.name)?'positive':'negative') as 'positive'|'negative',
        functionalGroup:atoms[0]!.residue.name}))]};
  };
  const protein = prepare("protein"), ligand = prepare("ligand");
  const waters = identities.flatMap(a => a?.entity === "water" && a.atomicNumber === 8 ? [{ oxygen: a }] : []);
  // Same metal list as native PLIP's supported metals / existing Web adapter.
  const symbols = new Map<number, string>([[3,"LI"],[11,"NA"],[12,"MG"],[13,"AL"],[19,"K"],[20,"CA"],[24,"CR"],[25,"MN"],[26,"FE"],[27,"CO"],[28,"NI"],[29,"CU"],[30,"ZN"],[31,"GA"],[37,"RB"],[38,"SR"],[44,"RU"],[45,"RH"],[46,"PD"],[47,"AG"],[48,"CD"],[49,"IN"],[51,"SB"],[55,"CS"],[56,"BA"],[57,"LA"],[58,"CE"],[59,"PR"],[62,"SM"],[63,"EU"],[64,"GD"],[65,"TB"],[70,"YB"],[71,"LU"],[74,"W"],[76,"OS"],[77,"IR"],[78,"PT"],[79,"AU"],[80,"HG"],[81,"TL"],[82,"PB"]]);
  return {
    protein,
    ligand: { ...ligand, metalTargets: [...ligand.metalTargets, ...waters.map(({ oxygen }) => ({ atom: oxygen, type: "O", location: "water" as const }))] },
    waters,
    metals: identities.flatMap(a => a && a.entity !== "water" && symbols.has(a.atomicNumber)
      ? [{ atom: a, element: symbols.get(a.atomicNumber)! }] : []),
  };
}
