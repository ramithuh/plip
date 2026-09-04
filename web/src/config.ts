/** PLIP 3.0.1 interaction thresholds. Distances are in Å and angles in degrees. */
export interface PlipThresholds {
  readonly minimumDistance: number;
  readonly hydrophobicDistanceMax: number;
  readonly hydrogenBondDistanceMax: number;
  readonly hydrogenBondDonorAngleMin: number;
  readonly piStackingDistanceMax: number;
  readonly piStackingAngleDeviation: number;
  readonly piStackingOffsetMax: number;
  readonly cationPiDistanceMax: number;
  readonly saltBridgeDistanceMax: number;
  readonly halogenDistanceMax: number;
  readonly halogenAcceptorAngle: number;
  readonly halogenDonorAngle: number;
  readonly halogenAngleDeviation: number;
  readonly waterBridgeDistanceMin: number;
  readonly waterBridgeDistanceMax: number;
  readonly waterBridgeOmegaMin: number;
  readonly waterBridgeOmegaMax: number;
  readonly waterBridgeThetaMin: number;
  readonly metalDistanceMax: number;
}

export const PLIP_DEFAULTS: PlipThresholds = Object.freeze({
  minimumDistance: 0.5,
  hydrophobicDistanceMax: 4.0,
  hydrogenBondDistanceMax: 4.1,
  hydrogenBondDonorAngleMin: 100,
  piStackingDistanceMax: 5.5,
  piStackingAngleDeviation: 30,
  piStackingOffsetMax: 2.0,
  cationPiDistanceMax: 6.0,
  saltBridgeDistanceMax: 5.5,
  halogenDistanceMax: 4.0,
  halogenAcceptorAngle: 120,
  halogenDonorAngle: 165,
  halogenAngleDeviation: 30,
  waterBridgeDistanceMin: 2.5,
  waterBridgeDistanceMax: 4.1,
  waterBridgeOmegaMin: 71,
  waterBridgeOmegaMax: 140,
  waterBridgeThetaMin: 100,
  metalDistanceMax: 3.0,
});
