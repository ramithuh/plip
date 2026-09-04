#!/usr/bin/env python3
"""Generate prepared-feature detector fixtures from this checkout's Python PLIP.

The script deliberately exercises PLIP's detection functions below the Open Babel
preparation boundary. Tiny Open Babel stubs provide only the atom/residue methods
read by those geometry functions; no chemical feature is inferred here.
"""

from __future__ import annotations

import json
import math
import hashlib
import sys
import types
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "python-plip-prepared.json"
TARGET_COMMIT = "2f4911d"


class PybelAtomStub:
    """Marker class used by plip.basic.supplemental.whichres*()."""


def _install_openbabel_import_stubs() -> None:
    openbabel_package = types.ModuleType("openbabel")
    openbabel_module = types.ModuleType("openbabel.openbabel")
    pybel_module = types.ModuleType("openbabel.pybel")
    openbabel_module.OBAtomAtomIter = lambda _atom: iter(())
    pybel_module.Atom = PybelAtomStub
    pybel_module.ob = openbabel_module
    openbabel_package.openbabel = openbabel_module
    openbabel_package.pybel = pybel_module
    sys.modules["openbabel"] = openbabel_package
    sys.modules["openbabel.openbabel"] = openbabel_module
    sys.modules["openbabel.pybel"] = pybel_module


_install_openbabel_import_stubs()
sys.path.insert(0, str(ROOT))

from plip.structure import detection  # noqa: E402


@dataclass(frozen=True)
class FakeResidue:
    name: str
    number: int
    chain: str

    def GetName(self) -> str:
        return self.name

    def GetNum(self) -> int:
        return self.number

    def GetChain(self) -> str:
        return self.chain

    def GetAtomProperty(self, _atom: Any, _property: int) -> bool:
        return True


class FakeOBAtom:
    def __init__(self, residue: FakeResidue, atom_type: str):
        self._residue = residue
        self._type = atom_type

    def GetResidue(self) -> FakeResidue:
        return self._residue

    def GetType(self) -> str:
        return self._type


class FakeAtom(PybelAtomStub):
    def __init__(self, payload: dict[str, Any]):
        residue = payload["residue"]
        self.idx = payload["index"]
        self.coords = tuple(payload["position"])
        self.atomicnum = payload["atomicNumber"]
        self.type = payload["name"]
        self.OBAtom = FakeOBAtom(
            FakeResidue(residue["name"], residue["number"], residue["chain"]),
            self.type,
        )


def _atom_cache(site: dict[str, Any]) -> dict[int, FakeAtom]:
    cache: dict[int, FakeAtom] = {}

    def visit(value: Any) -> None:
        if isinstance(value, dict):
            if {"index", "position", "residue", "entity"} <= value.keys():
                cache.setdefault(value["index"], FakeAtom(value))
            for child in value.values():
                visit(child)
        elif isinstance(value, list):
            for child in value:
                visit(child)

    visit(site)
    return cache


def _adapt(site: dict[str, Any]) -> SimpleNamespace:
    atoms = _atom_cache(site)

    def entity(payload: dict[str, Any]) -> SimpleNamespace:
        return SimpleNamespace(
            hydrophobic=[
                SimpleNamespace(atom=atoms[x["atom"]["index"]], orig_atom=atoms[x["atom"]["index"]], orig_idx=x["atom"]["parentIndex"])
                for x in payload["hydrophobic"]
            ],
            acceptors=[
                SimpleNamespace(a=atoms[x["atom"]["index"]], a_orig_atom=atoms[x["atom"]["index"]], a_orig_idx=x["atom"]["parentIndex"], type=x["type"])
                for x in payload["acceptors"]
            ],
            donors=[
                SimpleNamespace(
                    d=atoms[x["atom"]["index"]],
                    d_orig_atom=atoms[x["atom"]["index"]],
                    d_orig_idx=x["atom"]["parentIndex"],
                    h=atoms[x["hydrogen"]["index"]],
                    type=x["type"],
                )
                for x in payload["donors"]
            ],
            rings=[
                SimpleNamespace(
                    atoms=[atoms[a["index"]] for a in x["atoms"]],
                    orig_atoms=[atoms[a["index"]] for a in x["atoms"]],
                    atoms_orig_idx=[a["parentIndex"] for a in x["atoms"]],
                    normal=np.asarray(x["normal"], dtype=float),
                    center=np.asarray(x["center"], dtype=float),
                    obj=x["id"],
                )
                for x in payload["rings"]
            ],
            charges=[
                SimpleNamespace(
                    atoms=[atoms[a["index"]] for a in x["atoms"]],
                    orig_atoms=[atoms[a["index"]] for a in x["atoms"]],
                    atoms_orig_idx=[a["parentIndex"] for a in x["atoms"]],
                    center=np.asarray(x["center"], dtype=float),
                    type=x["charge"],
                    fgroup=x.get("functionalGroup"),
                    restype=x["atoms"][0]["residue"]["name"],
                    resnr=x["atoms"][0]["residue"]["number"],
                    reschain=x["atoms"][0]["residue"]["chain"],
                )
                for x in payload["charges"]
            ],
            halogen_acceptors=[
                SimpleNamespace(
                    o=atoms[x["oxygen"]["index"]],
                    o_orig_idx=x["oxygen"]["parentIndex"],
                    y=atoms[x["neighbor"]["index"]],
                    y_orig_idx=x["neighbor"]["parentIndex"],
                )
                for x in payload["halogenAcceptors"]
            ],
            halogen_donors=[
                SimpleNamespace(
                    x=atoms[x["halogen"]["index"]],
                    orig_x=atoms[x["halogen"]["index"]],
                    x_orig_idx=x["halogen"]["parentIndex"],
                    c=atoms[x["carbon"]["index"]],
                    c_orig_idx=[x["carbon"]["parentIndex"]],
                )
                for x in payload["halogenDonors"]
            ],
            metal_targets=[
                SimpleNamespace(
                    atom=atoms[x["atom"]["index"]],
                    atom_orig_idx=x["atom"]["parentIndex"],
                    type=x["type"],
                    restype=x["atom"]["residue"]["name"],
                    resnr=x["atom"]["residue"]["number"],
                    reschain=x["atom"]["residue"]["chain"],
                    location=x["location"],
                )
                for x in payload["metalTargets"]
            ],
        )

    return SimpleNamespace(
        protein=entity(site["protein"]),
        ligand=entity(site["ligand"]),
        waters=[SimpleNamespace(oxy=atoms[x["oxygen"]["index"]], oxy_orig_idx=x["oxygen"]["parentIndex"]) for x in site["waters"]],
        metals=[
            SimpleNamespace(m=atoms[x["atom"]["index"]], orig_m=atoms[x["atom"]["index"]], m_orig_idx=x["atom"]["parentIndex"])
            for x in site["metals"]
        ],
    )


def _residue(atom_payload: dict[str, Any]) -> dict[str, Any]:
    return atom_payload["residue"]


def _normalize(family: str, event: Any, protein_atom: dict[str, Any], ligand_atom: dict[str, Any], **extra: Any) -> dict[str, Any]:
    event_distance = getattr(event, "distance", None)
    if event_distance is None:
        event_distance = getattr(event, "distance_ad", None)
    if event_distance is None:
        event_distance = getattr(event, "distance_aw", None)
    if event_distance is None:
        raise ValueError(f"{family} oracle event has no supported distance field")
    output = {
        "family": family,
        "proteinResidue": _residue(protein_atom),
        "ligandResidue": _residue(ligand_atom),
        "distance": float(event_distance),
    }
    output.update(extra)
    return output


def _empty_entity() -> dict[str, list[Any]]:
    return {
        "hydrophobic": [], "acceptors": [], "donors": [], "rings": [], "charges": [],
        "halogenAcceptors": [], "halogenDonors": [], "metalTargets": [],
    }


_next_index = 1


def _atom(entity: str, position: list[float], *, name: str = "C", atomic_number: int = 6, residue_name: str | None = None) -> dict[str, Any]:
    global _next_index
    index = _next_index
    _next_index += 1
    protein = entity == "protein"
    water = entity == "water"
    return {
        "index": index,
        "parentIndex": index,
        "name": name,
        "atomicNumber": atomic_number,
        "position": position,
        "residue": {
            "name": residue_name or ("ALA" if protein else "HOH" if water else "LIG"),
            "number": 10 if protein else 100 if water else 1,
            "chain": "A" if protein else "W" if water else "L",
        },
        "entity": entity,
    }


def _site() -> dict[str, Any]:
    return {"protein": _empty_entity(), "ligand": _empty_entity(), "metals": [], "waters": []}


def _build_cases() -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []

    site = _site()
    p = _atom("protein", [0, 0, 0])
    l = _atom("ligand", [3, 0, 0])
    site["protein"]["hydrophobic"] = [{"atom": p, "neighbors": []}]
    site["ligand"]["hydrophobic"] = [{"atom": l, "neighbors": []}]
    adapted = _adapt(site)
    result = detection.hydrophobic_interactions(adapted.protein.hydrophobic, adapted.ligand.hydrophobic)
    cases.append({"name": "hydrophobic", "site": site, "expected": [_normalize("Hydrophobic", result[0], p, l)]})

    site = _site()
    p = _atom("protein", [2.8, 0, 0], name="O", atomic_number=8)
    d = _atom("ligand", [0, 0, 0], name="N", atomic_number=7)
    h = _atom("ligand", [1, 0, 0], name="H", atomic_number=1)
    site["protein"]["acceptors"] = [{"atom": p, "type": "O"}]
    site["ligand"]["donors"] = [{"atom": d, "hydrogen": h, "type": "N"}]
    adapted = _adapt(site)
    result = detection.hbonds(adapted.protein.acceptors, adapted.ligand.donors, False, "strong")
    cases.append({"name": "hydrogen-bond", "site": site, "expected": [_normalize("HydrogenBond", result[0], p, d, direction="ligand-donor")]})

    site = _site()
    p = _atom("protein", [0, 0, 0], residue_name="PHE")
    l = _atom("ligand", [0, 0, 4])
    site["protein"]["rings"] = [{"id": "p", "atoms": [p], "center": [0, 0, 0], "normal": [0, 0, 1]}]
    site["ligand"]["rings"] = [{"id": "l", "atoms": [l], "center": [0, 0, 4], "normal": [0.173648, 0, 0.984808]}]
    adapted = _adapt(site)
    result = detection.pistacking(adapted.protein.rings, adapted.ligand.rings)
    cases.append({"name": "pi-stacking", "site": site, "expected": [_normalize("PiStacking", result[0], p, l, subtype="parallel")]})

    site = _site()
    p = _atom("protein", [0, 0, 0], residue_name="TRP")
    l = _atom("ligand", [0, 0, 4], name="N", atomic_number=7)
    site["protein"]["rings"] = [{"id": "p", "atoms": [p], "center": [0, 0, 0], "normal": [0, 0, 1]}]
    site["ligand"]["charges"] = [{"atoms": [l], "center": [0, 0, 4], "charge": "positive"}]
    adapted = _adapt(site)
    result = detection.pication(adapted.protein.rings, adapted.ligand.charges, False)
    cases.append({"name": "cation-pi", "site": site, "expected": [_normalize("CationPi", result[0], p, l, direction="ligand-cation")]})

    site = _site()
    p = _atom("protein", [0, 0, 0], name="NZ", atomic_number=7, residue_name="LYS")
    l = _atom("ligand", [3, 0, 0], name="O", atomic_number=8)
    site["protein"]["charges"] = [{"atoms": [p], "center": [0, 0, 0], "charge": "positive"}]
    site["ligand"]["charges"] = [{"atoms": [l], "center": [3, 0, 0], "charge": "negative"}]
    adapted = _adapt(site)
    result = detection.saltbridge(adapted.protein.charges, adapted.ligand.charges, True)
    cases.append({"name": "salt-bridge", "site": site, "expected": [_normalize("SaltBridge", result[0], p, l, direction="protein-cation")]})

    site = _site()
    o = _atom("protein", [0, 0, 0], name="O", atomic_number=8)
    y = _atom("protein", [-0.5, math.sqrt(3) / 2, 0])
    x = _atom("ligand", [3, 0, 0], name="CL", atomic_number=17)
    c = _atom("ligand", [3 + math.cos(math.pi / 12), math.sin(math.pi / 12), 0])
    site["protein"]["halogenAcceptors"] = [{"oxygen": o, "neighbor": y, "type": "O"}]
    site["ligand"]["halogenDonors"] = [{"halogen": x, "carbon": c, "type": "halocarbon"}]
    adapted = _adapt(site)
    result = detection.halogen(adapted.protein.halogen_acceptors, adapted.ligand.halogen_donors)
    cases.append({"name": "halogen-bond", "site": site, "expected": [_normalize("HalogenBond", result[0], o, x)]})

    site = _site()
    water = _atom("water", [0, 0, 0], name="O", atomic_number=8)
    radians = math.radians(110)
    acceptor = _atom("ligand", [3 * math.cos(radians), 3 * math.sin(radians), 0], name="O", atomic_number=8)
    hydrogen = _atom("protein", [1, 0, 0], name="H", atomic_number=1)
    donor = _atom("protein", [3, 0, 0], name="N", atomic_number=7)
    site["ligand"]["acceptors"] = [{"atom": acceptor, "type": "O"}]
    site["protein"]["donors"] = [{"atom": donor, "hydrogen": hydrogen, "type": "N"}]
    site["waters"] = [{"oxygen": water}]
    adapted = _adapt(site)
    result = detection.water_bridges(
        adapted.protein.acceptors, adapted.ligand.acceptors,
        adapted.protein.donors, adapted.ligand.donors, adapted.waters,
    )
    cases.append({"name": "water-bridge", "site": site, "expected": [_normalize("WaterBridge", result[0], donor, acceptor, direction="protein-donor")]})

    site = _site()
    metal = _atom("ligand", [0, 0, 0], name="FE", atomic_number=26)
    target = _atom("protein", [2.2, 0, 0], name="NE2", atomic_number=7, residue_name="HIS")
    site["metals"] = [{"atom": metal, "element": "FE"}]
    site["protein"]["metalTargets"] = [{"atom": target, "type": "N", "location": "protein"}]
    adapted = _adapt(site)
    result = detection.metal_complexation(adapted.metals, adapted.ligand.metal_targets, adapted.protein.metal_targets)
    cases.append({"name": "metal-complex", "site": site, "expected": [_normalize("MetalComplex", result[0], target, metal, subtype="N")]})

    return cases


def main() -> None:
    detector_path = ROOT / "plip" / "structure" / "detection.py"
    document = {
        "oracle": {
            "implementation": "Python PLIP",
            "version": "3.0.1",
            "commit": TARGET_COMMIT,
            "detectorSha256": hashlib.sha256(detector_path.read_bytes()).hexdigest(),
        },
        "boundary": "prepared chemical features; Open Babel inference is deliberately excluded",
        "cases": _build_cases(),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    print(f"wrote {len(document['cases'])} Python-PLIP cases to {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
