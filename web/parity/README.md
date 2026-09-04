# Parity contract

PLIP Web has two distinct parity gates. Keeping them separate prevents a
correct distance calculation from being mistaken for a correct raw-structure
fingerprint.

## Gate A — prepared-feature detection

`generate_reference.py` calls the Python PLIP 3.0.1 detector functions in this
checkout and records their output for prepared hydrophobic atoms, hydrogen-bond
features, aromatic rings, charge groups, halogens, waters, and metals. The
TypeScript suite replays the same prepared inputs.

The checked-in fixture records the source commit and SHA-256 of
`plip/structure/detection.py`. Regenerate it with:

```bash
npm run parity:generate
npm test
```

Every supported prepared-feature case must match Python PLIP exactly on:

- interaction family;
- protein residue and chain;
- ligand residue and chain;
- donor/cation direction where applicable;
- PLIP subtype where applicable; and
- distance within floating-point tolerance.

Boundary and PLIP reporting-refinement behavior is covered separately by the
unit suite. This includes strict cutoffs, hydrophobic-patch reduction,
salt-bridge suppression of duplicate H-bonds, ring-contact suppression during
pi stacking, histidine cation-pi suppression, and water-bridge multiplicity.

## Gate B — raw PDB/mmCIF fingerprint

This gate is not passed yet. It requires a browser chemistry adapter that
matches the Open Babel preparation used by PLIP for ligand extraction,
connectivity, aromaticity, donor/acceptor assignment, charges, protonation, and
metal targets.

Before Weaver labels a result as PLIP, a deterministic raw-structure corpus
must be run through both Python PLIP and the browser adapter. Event comparison
uses the UI identity:

```text
focal ligand + normalized family + protein chain/residue + direction
```

Atom participation, distances, angles, and multiplicity remain secondary
fields and must also be reported when they disagree. Shared events sort first;
Python-only and browser-only events follow. Structures used for hydrogen-bond
parity must be protonated once and reused with Python PLIP's `--nohydro` mode so
Open Babel's hydrogen placement does not vary between runs.

Metal parity includes the selected coordination number, geometry label, RMS
score, and any target rejected by PLIP's geometry fit. This is important because
PLIP's greedy selection can produce a non-intuitive label even for an idealized
input; the browser port intentionally reproduces the reference behavior rather
than silently replacing it with a different classifier.

The raw gate should include positive and negative examples for every family,
multi-ligand and metal-containing systems, alternate locations, waters,
covalent ligands, nucleotides, and deliberately missing hydrogens. Until this
gate passes, the package compatibility metadata reports the chemistry backend
as caller-supplied.
