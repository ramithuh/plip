import { PLIP_DEFAULTS, type PlipThresholds } from "./config.js";
import { detectPreparedSite } from "./detection.js";
import type { DetectionResult, PreparedSite } from "./types.js";

/** Browser-side interaction detection over chemistry prepared by a caller. */
export class PlipFingerprint {
  constructor(readonly thresholds: PlipThresholds = PLIP_DEFAULTS) {}

  detect(site: PreparedSite): DetectionResult {
    return { events: detectPreparedSite(site, this.thresholds) };
  }
}
