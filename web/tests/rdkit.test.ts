import type { JSMol, RDKitLoader, RDKitModule } from "@rdkit/rdkit";
import { afterAll, beforeAll, describe, expect, it } from "vitest";

import {
  prepareRdkitEntity,
  prepareRdkitSite,
  type AtomRef,
  type RdkitComponentInput,
} from "../src/index.js";

interface RdkitJson {
  readonly defaults?: { readonly atom?: { readonly z?: number } };
  readonly molecules?: readonly { readonly atoms?: readonly { readonly z?: number }[] }[];
}

describe("RDKit.js feature preparation", () => {
  let rdkit: RDKitModule;
  let ownedMolecules: JSMol[];
  let nextAtomId: number;

  beforeAll(async () => {
    const namespace: unknown = await import("@rdkit/rdkit");
    const loader = (namespace as { readonly default?: RDKitLoader }).default;
    if (typeof loader !== "function") throw new TypeError("RDKit.js loader is unavailable");
    rdkit = await loader();
    ownedMolecules = [];
    nextAtomId = 1;
  });

  afterAll(() => ownedMolecules.forEach((molecule) => molecule.delete()));

  function component(
    smiles: string,
    entity: AtomRef["entity"] = "ligand",
    residueName = "LIG",
    atomNames: readonly string[] = [],
  ): RdkitComponentInput {
    const molecule = rdkit.get_mol(smiles, JSON.stringify({ removeHs: false, sanitize: true }));
    if (molecule === null) throw new Error(`RDKit could not parse ${smiles}`);
    ownedMolecules.push(molecule);
    const json = JSON.parse(molecule.get_json()) as RdkitJson;
    const rawAtoms = json.molecules?.[0]?.atoms ?? [];
    const defaultAtomicNumber = json.defaults?.atom?.z ?? 6;
    const atoms = rawAtoms.map((atom, index): AtomRef => {
      const id = nextAtomId++;
      const theta = 2 * Math.PI * index / Math.max(rawAtoms.length, 1);
      return {
        index: id,
        parentIndex: id,
        name: atomNames[index] ?? `A${index + 1}`,
        atomicNumber: atom.z ?? defaultAtomicNumber,
        position: [Math.cos(theta), Math.sin(theta), index % 2 === 0 ? 0 : 0.02],
        residue: { name: residueName, chain: entity === "protein" ? "A" : "L", number: 1 },
        entity,
      };
    });
    return { molecule, atoms };
  }

  it("prepares PLIP ligand charge, ring, halogen, and metal features", () => {
    const entities = [
      component("C[N+](C)(C)C"),
      component("CC(=O)[O-]"),
      component("CBr"),
      component("c1ccccc1"),
      component("c1ccccc1S"),
      component("[Mg+2]", "ligand", "MG"),
    ];
    const features = prepareRdkitEntity(rdkit, entities, "ligand");
    const site = prepareRdkitSite(rdkit, { protein: [], ligand: entities });

    expect(features.charges.map(({ functionalGroup }) => functionalGroup)).toEqual(
      expect.arrayContaining(["quartamine", "carboxylate"]),
    );
    expect(features.halogenDonors.map(({ type }) => type)).toContain("BR");
    expect(features.rings).toHaveLength(2);
    expect(features.metalTargets.map(({ type }) => type)).toContain("S");
    expect(site.metals.map(({ element }) => element)).toContain("MG");
  });

  it("does not invent donor geometry when only implicit hydrogens exist", () => {
    const implicit = prepareRdkitEntity(rdkit, [component("CN")], "ligand");
    const explicit = prepareRdkitEntity(rdkit, [component("[H]N([H])C")], "ligand");

    expect(implicit.donors).toHaveLength(0);
    expect(explicit.donors.length).toBeGreaterThan(0);
    expect(explicit.donors.every(({ hydrogen }) => hydrogen.atomicNumber === 1)).toBe(true);
  });

  it("uses protein residue identity for charge and metal-target typing", () => {
    const lysine = component("[NH4+]", "protein", "LYS", ["NZ"]);
    const aspartate = component("[O-]", "protein", "ASP", ["OD1"]);
    const features = prepareRdkitEntity(rdkit, [lysine, aspartate], "protein");

    expect(features.charges.map(({ charge, functionalGroup }) => `${functionalGroup}:${charge}`)).toEqual(
      expect.arrayContaining(["LYS:positive", "ASP:negative"]),
    );
    expect(features.metalTargets.map(({ atom }) => atom.name)).toContain("OD1");
  });

  it("only accepts oxygen atoms as water and exposes them as metal targets", () => {
    const waterO: AtomRef = {
      index: nextAtomId++, parentIndex: nextAtomId, name: "O", atomicNumber: 8,
      position: [0, 0, 0], residue: { name: "HOH", chain: "W", number: 1 }, entity: "water",
    };
    const waterH: AtomRef = {
      index: nextAtomId++, parentIndex: nextAtomId, name: "H1", atomicNumber: 1,
      position: [0.9, 0, 0], residue: { name: "HOH", chain: "W", number: 1 }, entity: "water",
    };
    const site = prepareRdkitSite(rdkit, { protein: [], ligand: [], waters: [waterO, waterH] });

    expect(site.waters.map(({ oxygen }) => oxygen.atomicNumber)).toEqual([8]);
    expect(site.ligand.metalTargets).toEqual([{ atom: waterO, type: "O", location: "water" }]);
  });
});
