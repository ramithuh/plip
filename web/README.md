# PLIP Web

`@ramithuh/plip-web` is a browser-compatible TypeScript port of PLIP's
interaction-detection stage. It targets PLIP 3.0.1 and retains PLIP's default
distance and angle thresholds.

This first milestone accepts a `PreparedSite`: explicit atoms plus the
hydrophobic, donor/acceptor, aromatic-ring, charge, halogen, water, and metal
features normally produced by PLIP's Open Babel preparation stage. It does not
yet claim end-to-end PLIP parity from a raw PDB or mmCIF file. In particular,
Open Babel chemical perception and protonation have not yet been ported to the
browser.

The separation is deliberate. Geometry can be tested independently while an
Open Babel/WASM or parity-validated replacement chemistry backend is developed.
Weaver should only expose an interaction family as “PLIP” after its raw-file
oracle tests pass.

The two-stage acceptance criteria and event-matching contract are documented in
[`parity/README.md`](parity/README.md).

## Development

```bash
npm install
npm run check
npm run test:coverage
```

The package is GPL-2.0-only, matching PLIP.
