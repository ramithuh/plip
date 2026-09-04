import type { Vec3 } from "./types.js";

export function distance(a: Vec3, b: Vec3): number {
  return Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

export function vector(from: Vec3, to: Vec3): Vec3 {
  return [to[0] - from[0], to[1] - from[1], to[2] - from[2]];
}

export function dot(a: Vec3, b: Vec3): number {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

export function norm(a: Vec3): number {
  return Math.hypot(a[0], a[1], a[2]);
}

export function angle(a: Vec3, b: Vec3): number {
  if (a[0] === b[0] && a[1] === b[1] && a[2] === b[2]) return 0;
  const denominator = norm(a) * norm(b);
  if (denominator === 0) return Number.NaN;
  const cosine = Math.max(-1, Math.min(1, dot(a, b) / denominator));
  return Math.acos(cosine) * (180 / Math.PI);
}

/** Orthogonal projection of `point` onto the plane through `planePoint`. */
export function projectOntoPlane(normal: Vec3, planePoint: Vec3, point: Vec3): Vec3 {
  const denominator = dot(normal, normal);
  if (denominator === 0) return point;
  const scale = -dot(normal, vector(planePoint, point)) / denominator;
  return [
    point[0] + scale * normal[0],
    point[1] + scale * normal[1],
    point[2] + scale * normal[2],
  ];
}
