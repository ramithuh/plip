import { describe, expect, it } from "vitest";
import { angle, distance, dot, norm, projectOntoPlane, vector } from "../src/index.js";

describe("geometry primitives", () => {
  it("matches elementary Euclidean vector operations", () => {
    expect(distance([0, 0, 0], [2, 3, 6])).toBe(7);
    expect(vector([1, 2, 3], [4, 6, 3])).toEqual([3, 4, 0]);
    expect(dot([1, 2, 3], [4, 5, 6])).toBe(32);
    expect(norm([3, 4, 0])).toBe(5);
  });

  it("returns stable angles and handles degenerate vectors", () => {
    expect(angle([1, 0, 0], [0, 1, 0])).toBeCloseTo(90);
    expect(angle([1, 0, 0], [-1, 0, 0])).toBeCloseTo(180);
    expect(angle([1, 0, 0], [1, 0, 0])).toBe(0);
    expect(angle([0, 0, 0], [1, 0, 0])).toBeNaN();
  });

  it("projects a point onto an arbitrary plane", () => {
    expect(projectOntoPlane([0, 0, 2], [1, 1, 1], [4, 5, 8])).toEqual([4, 5, 1]);
    expect(projectOntoPlane([0, 0, 0], [1, 1, 1], [4, 5, 8])).toEqual([4, 5, 8]);
  });
});
