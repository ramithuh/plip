import { describe, expect, it } from "vitest";
import reference from "./fixtures/python-plip-prepared.json";
import { detectPreparedSite, type InteractionEvent, type PreparedSite } from "../src/index.js";

interface ExpectedEvent {
  readonly family: InteractionEvent["family"];
  readonly proteinResidue: InteractionEvent["proteinResidue"];
  readonly ligandResidue: InteractionEvent["ligandResidue"];
  readonly distance: number;
  readonly direction?: InteractionEvent["direction"];
  readonly subtype?: string;
  readonly geometry?: Readonly<Record<string, number | string>>;
}

function normalize(event: InteractionEvent): ExpectedEvent {
  return {
    family: event.family,
    proteinResidue: event.proteinResidue,
    ligandResidue: event.ligandResidue,
    distance: event.distance,
    ...(event.direction === undefined ? {} : { direction: event.direction }),
    ...(event.subtype === undefined ? {} : { subtype: event.subtype }),
    ...(event.geometry === undefined ? {} : { geometry: event.geometry }),
  };
}

describe("Python PLIP prepared-feature oracle", () => {
  it("is pinned to the PLIP version targeted by this port", () => {
    expect(reference.oracle.version).toBe("3.0.1");
    expect(reference.boundary).toContain("Open Babel inference is deliberately excluded");
  });

  for (const fixture of reference.cases) {
    it(`matches Python PLIP for ${fixture.name}`, () => {
      const actual = detectPreparedSite(fixture.site as PreparedSite).map(normalize);
      expect(actual).toHaveLength(fixture.expected.length);
      for (const [index, expected] of fixture.expected.entries()) {
        expect(actual[index]).toMatchObject({
          family: expected.family,
          proteinResidue: expected.proteinResidue,
          ligandResidue: expected.ligandResidue,
          ...(expected.direction === undefined ? {} : { direction: expected.direction }),
          ...(expected.subtype === undefined ? {} : { subtype: expected.subtype }),
          ...(expected.geometry === undefined ? {} : {
            geometry: {
              coordination: expected.geometry.coordination,
              observedCoordination: expected.geometry.observedCoordination,
              shape: expected.geometry.shape,
            },
          }),
        });
        expect(actual[index]!.distance).toBeCloseTo(expected.distance, 6);
        if (expected.geometry?.rms !== undefined) {
          expect(actual[index]!.geometry?.rms).toBeCloseTo(expected.geometry.rms as number, 6);
        }
      }
    });
  }
});
