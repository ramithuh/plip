# PLIP Web

`@ramithuh/plip-web` is a browser-compatible TypeScript port of PLIP's
interaction-detection stage. It targets PLIP 3.0.1 and retains PLIP's default
distance and angle thresholds.

The exact detector accepts a `PreparedSite`: explicit atoms plus the
hydrophobic, donor/acceptor, aromatic-ring, charge, halogen, water, and metal
features normally produced by PLIP's Open Babel preparation stage. An
experimental `prepareRdkitSite()` adapter now derives those features from
sanitized RDKit.js components while retaining the caller's atom identities and
3D coordinates.

`prepareOpenBabelSite(molecule, identities)` now also accepts the JSON graph and
perception flags from original Open Babel C++ compiled to WASM. It uses native
donor/acceptor flags, hybridization, implicit-H counts, bonds and rings, while
sharing the Web port's existing functional-group and interaction feature rules.
It does not call RDKit to retype the Open Babel output. Entries in `identities`
may be undefined for contextual atoms outside the selected site.

Weaver now uses this backend by default, prepares the complete parsed structure
in a worker, and offers polar-H addition or input-H-only modes independently of
ProLIF. The WASM bridge/build sources and runtime assets live in Weaver's
`vendor/openbabel-wasm` and `public/openbabel`. Atom indices returned by the
bridge are 1-based OB indices; the caller must validate and map them back to
original structure identities (including generated H) before using this API.

This is still fixture-bounded feature adaptation, not a full port of Python
PLIP's ligand selection and preprocessing. Weaver's fixed-protonated 4DST test
matches all six native PLIP 3.0.1 atom-level events. Open Babel can generate
variable H orientations for underdetermined geometry even in native runs;
independently protonated results need a separate sensitivity comparison.

The package deliberately does not own PDB/mmCIF parsing or receptor/ligand
selection. A caller supplies one RDKit molecule and ordered `AtomRef` array per
component. Their atom counts, atomic numbers, and order are validated before
features are prepared. This makes the chemistry adapter reusable without
silently changing structure semantics.

The RDKit adapter is not yet a universal substitute for PLIP's Open Babel
preparation. It requires explicit donor hydrogens for hydrogen-bond geometry,
and its end-to-end parity evidence is currently a bounded corpus. Weaver should
label results with that measured compatibility rather than imply all raw
structures have passed.

The prepared-feature core covers PLIP's hydrophobic contacts, hydrogen bonds,
pi stacking, cation-pi contacts, salt bridges, halogen bonds, water bridges,
and metal complexes. Metal targets are fitted with PLIP's original greedy
coordination-number and geometry heuristic before cross-interface contacts are
reported.

The two-stage acceptance criteria and event-matching contract are documented in
[`parity/README.md`](parity/README.md).

## RDKit.js adapter

```ts
import initRDKitModule from "@rdkit/rdkit";
import { detectPreparedSite, prepareRdkitSite } from "@ramithuh/plip-web";

const rdkit = await initRDKitModule();
const site = prepareRdkitSite(rdkit, {
  protein: [{ molecule: proteinMol, atoms: proteinAtoms }],
  ligand: [{ molecule: ligandMol, atoms: ligandAtoms }],
  waters,
});
const events = detectPreparedSite(site);
```

`molecule` must be sanitized, must preserve explicit hydrogens, and must have
the same atom order as `atoms`. The original coordinates in `atoms` are used for
every reported distance and angle.

## Development

```bash
npm install
npm run check
npm run test:coverage
```

The package is GPL-2.0-only, matching PLIP.
