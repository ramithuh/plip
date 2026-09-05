import { describe, expect, it } from "vitest";
import { prepareOpenBabelSite, type AtomRef, type OpenBabelMolecule } from "../src/index.js";

const molecule: OpenBabelMolecule = {
  ok: true, version: "3.2.1", atomsBefore: 2, atomsAfter: 3, addedHydrogens: 1,
  atoms: [
    { index: 1, element: 7, formalCharge: 0, hybridization: 3, implicitHydrogens: 0, aromatic: false, donor: true, acceptor: false, donorHydrogen: false, xyz: [0,0,0] },
    { index: 2, element: 8, formalCharge: 0, hybridization: 2, implicitHydrogens: 0, aromatic: false, donor: false, acceptor: true, donorHydrogen: false, xyz: [3,0,0] },
    { index: 3, element: 1, formalCharge: 0, hybridization: 0, implicitHydrogens: 0, aromatic: false, donor: false, acceptor: false, donorHydrogen: true, xyz: [1,0,0] },
  ],
  bonds: [{ begin: 1, end: 3, order: 1, aromatic: false }], rings: [],
};
const refs: AtomRef[] = molecule.atoms.map(a => ({ index: a.index, parentIndex: a.index,
  atomicNumber: a.element, position: a.xyz, name: String(a.index),
  entity: a.index === 2 ? "protein" : "ligand", residue: { name: a.index === 2 ? "SER" : "LIG", chain: "A", number: a.index === 2 ? 1 : 2 } }));

describe("native Open Babel feature adapter", () => {
  it('uses native residue-property charge membership, not atom-name guesses', () => {
    const m: OpenBabelMolecule = {...molecule, atomsBefore:2,atomsAfter:2,addedHydrogens:0,bonds:[],
      atoms:molecule.atoms.slice(0,2).map(a=>({...a,element:7,type:'N3',residueProperty8:true}))};
    const ids: AtomRef[] = m.atoms.map((a,i)=>({...refs[i]!,atomicNumber:7,entity:'protein',
      name:i===0?'N':'NZ',residue:{name:'LYS',chain:'A',number:1}}));
    expect(prepareOpenBabelSite(m,ids).protein.charges[0]!.atoms.map(a=>a.name)).toEqual(['N','NZ']);
  });
  it('classifies ligand charge from heavy neighbors, separately from donor H',()=>{
    const m: OpenBabelMolecule={...molecule,atomsBefore:4,atomsAfter:4,addedHydrogens:0,rings:[],
      atoms:[7,6,6,1].map((element,i)=>({...molecule.atoms[0]!,index:i+1,element})),
      bonds:[2,3,4].map(end=>({begin:1,end,order:1,aromatic:false}))};
    const ids: AtomRef[]=m.atoms.map(a=>({...refs[0]!,index:a.index,parentIndex:a.index,atomicNumber:a.element}));
    expect(prepareOpenBabelSite(m,ids).ligand.charges).toEqual([]);
  });
  it('orders ring probe atoms by original index as native PLIP does',()=>{
    const m: OpenBabelMolecule={...molecule,atomsBefore:6,atomsAfter:6,addedHydrogens:0,
      atoms:Array.from({length:6},(_,i)=>({...molecule.atoms[0]!,index:i+1,element:6,aromatic:true,
        xyz:[Math.cos(i*Math.PI/3),Math.sin(i*Math.PI/3),0] as const})),bonds:[],
      rings:[{atoms:[4,3,2,1,6,5],aromatic:true}]};
    const ids: AtomRef[]=m.atoms.map(a=>({...refs[0]!,index:a.index,parentIndex:a.index,atomicNumber:6,position:a.xyz}));
    expect(prepareOpenBabelSite(m,ids).ligand.rings[0]!.atoms.map(a=>a.parentIndex)).toEqual([1,2,3,4,5,6]);
  });
  it("uses native donor/acceptor flags and native donor-H identity", () => {
    const site = prepareOpenBabelSite(molecule, refs);
    expect(site.ligand.acceptors).toHaveLength(0);
    expect(site.protein.acceptors[0]!.atom).toBe(refs[1]);
    expect(site.ligand.donors[0]!.hydrogen).toBe(refs[2]);
    const noDonorH = { ...molecule, atoms: molecule.atoms.map(a => ({ ...a, donorHydrogen: false })) };
    expect(prepareOpenBabelSite(noDonorH, refs).ligand.donors).toHaveLength(0);
  });
  it("does not treat unselected contextual atoms as receptor or ligand", () => {
    const site = prepareOpenBabelSite(molecule, [undefined, refs[1], undefined]);
    expect(site.ligand.donors).toHaveLength(0);
    expect(site.protein.acceptors).toHaveLength(1);
  });
  it("rejects count and element mismatches", () => {
    expect(() => prepareOpenBabelSite(molecule, [])).toThrow(/count/);
    expect(() => prepareOpenBabelSite(molecule, [{ ...refs[0]!, atomicNumber: 6 }, ...refs.slice(1)])).toThrow(/element/);
  });
  it("rejects malformed bond endpoints", () => {
    expect(() => prepareOpenBabelSite({ ...molecule, bonds: [{ begin: 1, end: 99, order: 1, aromatic: false }] }, refs)).toThrow(/endpoint/);
  });
});
