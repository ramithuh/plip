export const PLIP_WEB_COMPATIBILITY = Object.freeze({
  targetVersion: "3.0.1",
  targetCommit: "2f4911d",
  status: "prepared detection parity; experimental raw-input parity",
  chemistryBackend: "RDKit.js adapter over caller-parsed components",
  limitations: [
    "PDB/mmCIF parsing, component selection, and atom-identity mapping remain caller responsibilities.",
    "RDKit.js does not universally reproduce PLIP's Open Babel chemical perception; raw parity is currently fixture-bounded.",
    "Hydrogen-bond and water-bridge detection requires explicit donor hydrogens in the input structure.",
    "Metal contact detection and PLIP's coordination-geometry fitting are implemented over prepared targets.",
  ],
});
