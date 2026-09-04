import type { JSMol, RDKitLoader, RDKitModule } from "@rdkit/rdkit";
import { createHash } from "node:crypto";
import { readFileSync } from "node:fs";
import { afterAll, beforeAll, describe, expect, it } from "vitest";

import {
  detectPreparedSite,
  prepareRdkitSite,
  type AtomRef,
  type InteractionEvent,
  type RdkitComponentInput,
} from "../src/index.js";

const ATOMIC_NUMBER = new Map([
  ["H", 1], ["C", 6], ["N", 7], ["O", 8], ["F", 9], ["MG", 12],
  ["P", 15], ["S", 16], ["CL", 17], ["BR", 35], ["I", 53],
]);
const ELEMENT_SYMBOL = new Map([...ATOMIC_NUMBER].map(([symbol, number]) => [number, symbol]));
const SOURCE = readFileSync(
  new URL("../../plip/test/pdb/4dst_protonated.pdb", import.meta.url),
  "utf8",
);
const LINES = SOURCE.split(/\r?\n/);

interface SerializedComponent {
  readonly molblock: string;
  readonly atoms: readonly AtomRef[];
}

function atomFromPdb(line: string, entity: AtomRef["entity"]): AtomRef {
  const rawElement = line.slice(76, 78).trim()
    || line.slice(12, 16).trim().replace(/[0-9]/g, "").slice(0, 2);
  const element = rawElement.toUpperCase();
  const atomicNumber = ATOMIC_NUMBER.get(element);
  if (atomicNumber === undefined) throw new Error(`Unsupported test element ${element}`);
  const serial = Number.parseInt(line.slice(6, 11), 10);
  return {
    index: serial,
    parentIndex: serial,
    name: line.slice(12, 16).trim(),
    atomicNumber,
    position: [
      Number.parseFloat(line.slice(30, 38)),
      Number.parseFloat(line.slice(38, 46)),
      Number.parseFloat(line.slice(46, 54)),
    ],
    residue: {
      name: line.slice(17, 20).trim(),
      chain: line.slice(21, 22).trim(),
      number: Number.parseInt(line.slice(22, 26), 10),
      insertionCode: line.slice(26, 27).trim(),
    },
    entity,
  };
}

function conectMultiplicity(selectedSerials: ReadonlySet<number>): Map<string, number> {
  const directed = new Map<string, number>();
  for (const line of LINES) {
    if (!line.startsWith("CONECT")) continue;
    const serials = line.slice(6).match(/.{1,5}/g)
      ?.map((chunk) => Number.parseInt(chunk.trim(), 10))
      .filter(Number.isFinite) ?? [];
    const source = serials[0];
    if (source === undefined || !selectedSerials.has(source)) continue;
    for (const target of serials.slice(1)) {
      if (!selectedSerials.has(target)) continue;
      const key = `${source}:${target}`;
      directed.set(key, (directed.get(key) ?? 0) + 1);
    }
  }
  return directed;
}

function molblockFor(
  predicate: (line: string) => boolean,
  entity: AtomRef["entity"],
): SerializedComponent {
  const atomLines = LINES.filter((line) =>
    (line.startsWith("ATOM  ") || line.startsWith("HETATM")) && predicate(line));
  const atoms = atomLines.map((line) => atomFromPdb(line, entity));
  const selected = new Set(atoms.map((atom) => atom.index));
  const localIndex = new Map(atoms.map((atom, index) => [atom.index, index + 1]));
  const directed = conectMultiplicity(selected);
  const bonds: Array<readonly [number, number, number]> = [];
  for (const source of selected) {
    for (const target of selected) {
      if (source >= target) continue;
      const order = Math.max(
        directed.get(`${source}:${target}`) ?? 0,
        directed.get(`${target}:${source}`) ?? 0,
      );
      if (order > 0) bonds.push([
        localIndex.get(source)!, localIndex.get(target)!, Math.min(order, 3),
      ]);
    }
  }
  const atomBlock = atoms.map((atom) => {
    const symbol = ELEMENT_SYMBOL.get(atom.atomicNumber)!;
    const normalizedSymbol = symbol[0]! + symbol.slice(1).toLowerCase();
    return [
      atom.position[0].toFixed(4).padStart(10),
      atom.position[1].toFixed(4).padStart(10),
      atom.position[2].toFixed(4).padStart(10),
      ` ${normalizedSymbol.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`,
    ].join("");
  });
  const bondBlock = bonds.map(([left, right, order]) =>
    `${String(left).padStart(3)}${String(right).padStart(3)}${String(order).padStart(3)}  0  0  0  0`);
  return {
    atoms,
    molblock: [
      "", "  PLIP Web raw-parity fixture", "",
      `${String(atoms.length).padStart(3)}${String(bonds.length).padStart(3)}  0  0  0  0            999 V2000`,
      ...atomBlock, ...bondBlock, "M  END",
    ].join("\n"),
  };
}

function rdkitComponent(rdkit: RDKitModule, serialized: SerializedComponent): RdkitComponentInput {
  const molecule = rdkit.get_mol(serialized.molblock, JSON.stringify({
    removeHs: false,
    sanitize: true,
  }));
  if (molecule === null) throw new Error("RDKit could not parse a raw-parity component");
  return { molecule, atoms: serialized.atoms };
}

interface EventIdentity {
  readonly family: InteractionEvent["family"];
  readonly protein: string;
  readonly ligand: string;
  readonly direction: string;
  readonly proteinAtoms: readonly number[];
  readonly ligandAtoms: readonly number[];
  readonly water?: number;
}

const REFERENCE = JSON.parse(readFileSync(
  new URL("./fixtures/python-plip-raw-4dst.json", import.meta.url),
  "utf8",
)) as {
  readonly oracle: { readonly input_sha256: string; readonly site: string };
  readonly events: readonly EventIdentity[];
};

function eventIdentity(event: InteractionEvent): EventIdentity {
  const residue = (value: InteractionEvent["proteinResidue"]): string =>
    `${value.name}:${value.chain}:${value.number}${value.insertionCode ?? ""}`;
  return {
    family: event.family,
    protein: residue(event.proteinResidue),
    ligand: residue(event.ligandResidue),
    direction: event.direction ?? "",
    proteinAtoms: event.proteinAtoms.map((atom) => atom.parentIndex).sort((a, b) => a - b),
    ligandAtoms: event.ligandAtoms.map((atom) => atom.parentIndex).sort((a, b) => a - b),
    ...(event.water === undefined ? {} : { water: event.water.parentIndex }),
  };
}

function identitySort(left: EventIdentity, right: EventIdentity): number {
  return JSON.stringify(left).localeCompare(JSON.stringify(right));
}

describe("RDKit.js raw-structure preparation parity", () => {
  let rdkit: RDKitModule;
  let ownedMolecules: JSMol[];

  beforeAll(async () => {
    const namespace: unknown = await import("@rdkit/rdkit");
    const loader = (namespace as { readonly default?: RDKitLoader }).default;
    if (typeof loader !== "function") throw new TypeError("RDKit.js loader is unavailable");
    rdkit = await loader();
    ownedMolecules = [];
  });

  afterAll(() => {
    ownedMolecules.forEach((molecule) => molecule.delete());
  });

  it("matches Python PLIP 3.0.1 on its protonated 4DST/9LI fixture", () => {
    expect(createHash("sha256").update(SOURCE).digest("hex")).toBe(REFERENCE.oracle.input_sha256);
    expect(REFERENCE.oracle.site).toBe("9LI:A:201");
    const proteinAtomLines = LINES.filter((line) => line.startsWith("ATOM  "));
    const residueKeys = [...new Set(proteinAtomLines.map((line) =>
      `${line.slice(21, 22)}:${line.slice(22, 27)}`))];
    const protein = residueKeys.map((key) => rdkitComponent(rdkit, molblockFor(
      (line) => line.startsWith("ATOM  ")
        && `${line.slice(21, 22)}:${line.slice(22, 27)}` === key,
      "protein",
    )));
    const ligand = rdkitComponent(rdkit, molblockFor(
      (line) => line.slice(17, 20).trim() === "9LI"
        && line.slice(21, 22) === "A"
        && Number.parseInt(line.slice(22, 26), 10) === 201,
      "ligand",
    ));
    ownedMolecules.push(...protein.map(({ molecule }) => molecule), ligand.molecule);
    const waters = LINES.filter((line) =>
      line.startsWith("HETATM") && line.slice(17, 20).trim() === "HOH")
      .map((line) => atomFromPdb(line, "water"));

    const actual = detectPreparedSite(prepareRdkitSite(rdkit, {
      protein,
      ligand: [ligand],
      waters,
    })).map(eventIdentity).sort(identitySort);
    const expected = [...REFERENCE.events].sort(identitySort);

    expect(actual).toEqual(expected);
    expect(new Set(actual.map((event) => event.family))).toEqual(new Set([
      "Hydrophobic", "HydrogenBond", "CationPi", "WaterBridge",
    ]));
  });

  it("rejects a structure identity array that does not follow RDKit atom order", () => {
    const serialized = molblockFor(
      (line) => line.slice(17, 20).trim() === "9LI"
        && line.slice(21, 22) === "A"
        && Number.parseInt(line.slice(22, 26), 10) === 201,
      "ligand",
    );
    const component = rdkitComponent(rdkit, serialized);
    ownedMolecules.push(component.molecule);
    const wrongAtoms = component.atoms.map((atom, index) =>
      index === 0 ? { ...atom, atomicNumber: 7 } : atom);
    expect(() => prepareRdkitSite(rdkit, {
      protein: [],
      ligand: [{ molecule: component.molecule, atoms: wrongAtoms }],
    })).toThrow("supplied identity");
  });
});
