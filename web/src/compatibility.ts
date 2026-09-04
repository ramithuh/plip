export const PLIP_WEB_COMPATIBILITY = Object.freeze({
  targetVersion: "3.0.1",
  targetCommit: "2f4911d",
  status: "prepared-feature detection core",
  chemistryBackend: "caller-supplied",
  limitations: [
    "Open Babel structure preparation and chemical perception are not yet browser-ported.",
    "Hydrogen-bond and water-bridge detection requires explicit prepared donor hydrogens.",
    "Metal contact detection and PLIP's coordination-geometry fitting are implemented over prepared targets.",
  ],
});
