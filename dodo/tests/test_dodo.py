"""
Unit and regression test for the dodo package.
"""

# Import package, test suite, and other packages as needed
import sys
from types import SimpleNamespace

import pytest

import dodo
from dodo.pdb_tools import PDBParser
from dodo.pulchra_tools import rebuild_PDBParserObj_with_pulchra


def _write_fake_pulchra(tmp_path, atoms):
    fake_pulchra = tmp_path / "pulchra"
    fake_pulchra.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "from pathlib import Path\n"
        f"ATOMS = {atoms!r}\n"
        "def atom_line(serial, atom, residue, residue_index, x, y, z):\n"
        "    return f'ATOM  {serial:5d} {atom:>4s} {residue:<3s} A{residue_index:4d}    {x:8.3f}{y:8.3f}{z:8.3f}      {1.0:6.2f}\\n'\n"
        "input_pdb = Path(sys.argv[-1])\n"
        "rebuilt_pdb = input_pdb.with_name(f'{input_pdb.stem}.rebuilt.pdb')\n"
        "rebuilt_pdb.write_text(''.join(atom_line(*atom) for atom in ATOMS))\n"
    )
    fake_pulchra.chmod(0o755)
    return fake_pulchra


def test_dodo_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "dodo" in sys.modules


def test_pulchra_wrapper_uses_rebuilt_output(tmp_path):
    fake_pulchra = _write_fake_pulchra(tmp_path, [
        (1, "N", "ALA", 1, 9.0, 9.0, 9.0),
        (2, "CA", "ALA", 1, 9.1, 9.1, 9.1),
        (3, "C", "ALA", 1, 9.2, 9.2, 9.2),
        (4, "N", "GLY", 2, 5.0, 5.0, 5.0),
        (5, "CA", "GLY", 2, 5.1, 5.1, 5.1),
        (6, "C", "GLY", 2, 5.2, 5.2, 5.2),
    ])

    pdb_obj = SimpleNamespace(
        all_atom_coords_by_index={
            0: {"N": (1.0, 1.0, 1.0), "CA": (2.0, 2.0, 2.0), "C": (3.0, 3.0, 3.0)},
            1: {"CA": (3.8, 0.0, 0.0)},
        },
        all_atom_coords_by_index_with_aa={
            0: {"amino_acid": "ALA", "coords": {"N": (1.0, 1.0, 1.0), "CA": (2.0, 2.0, 2.0), "C": (3.0, 3.0, 3.0)}},
            1: {"amino_acid": "GLY", "coords": {"CA": (3.8, 0.0, 0.0)}},
        },
        index_to_3aa={0: "ALA", 1: "GLY"},
        sequence="AG",
        sequence_3aa_by_index={0: {"N": "ALA", "CA": "ALA", "C": "ALA"}, 1: {"CA": "GLY"}},
        beta_vals_by_index={0: 1.0, 1: 1.0},
        regions_dict={"folded_1": [0, 0], "idr_2": [1, 1]},
        number_atoms=4,
    )
    out_path = tmp_path / "rebuilt.pdb"

    rebuild_PDBParserObj_with_pulchra(pdb_obj, out_path, executable=fake_pulchra)
    rebuilt = PDBParser(out_path.read_text().splitlines())
    ca_elements = [
        line[76:78].strip()
        for line in out_path.read_text().splitlines()
        if line.startswith("ATOM  ") and line[12:16].strip() == "CA"
    ]

    assert rebuilt.all_atom_coords_by_index[-1]["N"] == (1.0, 1.0, 1.0)
    assert rebuilt.all_atom_coords_by_index[0]["N"] == (5.0, 5.0, 5.0)
    assert ca_elements == ["C", "C"]


def test_pulchra_wrapper_maps_trace_indices_to_original_indices(tmp_path):
    fake_pulchra = _write_fake_pulchra(tmp_path, [
        (1, "N", "MET", 1, 8.0, 8.0, 8.0),
        (2, "CA", "MET", 1, 8.1, 8.1, 8.1),
        (3, "C", "MET", 1, 8.2, 8.2, 8.2),
        (4, "N", "GLY", 2, 5.0, 5.0, 5.0),
        (5, "CA", "GLY", 2, 5.1, 5.1, 5.1),
        (6, "C", "GLY", 2, 5.2, 5.2, 5.2),
    ])

    pdb_obj = SimpleNamespace(
        all_atom_coords_by_index={
            -1: {"N": (1.0, 1.0, 1.0), "CA": (2.0, 2.0, 2.0), "C": (3.0, 3.0, 3.0)},
            0: {"CA": (3.8, 0.0, 0.0)},
        },
        all_atom_coords_by_index_with_aa={
            -1: {"amino_acid": "MET", "coords": {"N": (1.0, 1.0, 1.0), "CA": (2.0, 2.0, 2.0), "C": (3.0, 3.0, 3.0)}},
            0: {"amino_acid": "GLY", "coords": {"CA": (3.8, 0.0, 0.0)}},
        },
        index_to_3aa={-1: "MET", 0: "GLY"},
        sequence="MG",
        sequence_3aa_by_index={-1: {"N": "MET", "CA": "MET", "C": "MET"}, 0: {"CA": "GLY"}},
        beta_vals_by_index={-1: 1.0, 0: 1.0},
        regions_dict={"folded_1": [-1, -1], "idr_2": [0, 0]},
        number_atoms=4,
    )
    out_path = tmp_path / "rebuilt.pdb"

    rebuild_PDBParserObj_with_pulchra(pdb_obj, out_path, executable=fake_pulchra)
    gly_n_line = next(
        line for line in out_path.read_text().splitlines()
        if line.startswith("ATOM  ")
        and line[12:16].strip() == "N"
        and line[17:20].strip() == "GLY"
    )

    assert gly_n_line[22:26].strip() == "0"
    assert (
        float(gly_n_line[30:38]),
        float(gly_n_line[38:46]),
        float(gly_n_line[46:54]),
    ) == (5.0, 5.0, 5.0)
