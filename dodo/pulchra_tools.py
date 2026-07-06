import shutil
import subprocess
import tempfile
from pathlib import Path

from dodo.dodo_exceptions import dodoException
from dodo.pdb_tools import PDBParser, save_pdb_from_PDBParserObj


def _resolve_pulchra_executable(executable):
    executable_path = Path(executable).expanduser()
    if executable_path.parent != Path(".") or executable_path.exists():
        if not executable_path.is_file():
            raise dodoException(f'PULCHRA executable does not exist: {executable}')
        return str(executable_path.resolve())

    resolved = shutil.which(executable)
    if resolved is None:
        raise dodoException(
            f'Could not find PULCHRA executable "{executable}". '
            'Install PULCHRA or pass pulchra_executable with the full path.'
        )
    return resolved


def _rebuilt_path(input_pdb):
    input_pdb = Path(input_pdb)
    return input_pdb.with_name(f'{input_pdb.stem}.rebuilt.pdb')


def _protein_atom_element(atom_name):
    atom_name = atom_name.strip()
    if atom_name == '':
        return '  '
    if atom_name[0].isdigit():
        atom_name = atom_name[1:]
    return atom_name[0].upper().rjust(2)


def _normalize_protein_atom_elements(pdb_path):
    pdb_path = Path(pdb_path)
    fixed_lines = []
    for line in pdb_path.read_text().splitlines(keepends=True):
        if not line.startswith('ATOM  '):
            fixed_lines.append(line)
            continue

        newline = '\n' if line.endswith('\n') else ''
        atom_line = line.rstrip('\n').ljust(78)
        element = _protein_atom_element(atom_line[12:16])
        fixed_lines.append(f'{atom_line[:76]}{element}{atom_line[78:]}{newline}')

    pdb_path.write_text(''.join(fixed_lines))


def _format_ca_line(atom_index, residue_name, residue_index, xyz, beta):
    x, y, z = xyz
    return (
        f'ATOM  {atom_index:5d}  CA  {residue_name:>3s} A{residue_index:4d}    '
        f'{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{float(beta):6.2f}           C  \n'
    )


def _write_ca_trace(ca_records, out_path):
    Path(out_path).write_text(''.join(
        _format_ca_line(
            atom_index=i + 1,
            residue_name=record['residue_name'],
            residue_index=i + 1,
            xyz=record['xyz'],
            beta=record['beta'],
        )
        for i, record in enumerate(ca_records)
    ))


def run_pulchra(input_pdb, out_path=None, executable='pulchra', verbose=False):
    input_pdb = Path(input_pdb)
    pulchra_executable = _resolve_pulchra_executable(executable)
    rebuilt_path = _rebuilt_path(input_pdb)

    if rebuilt_path.exists():
        rebuilt_path.unlink()

    command = [pulchra_executable, input_pdb.name]
    if verbose:
        print(f'Running PULCHRA: {" ".join(command)}')

    result = subprocess.run(command, cwd=str(input_pdb.parent), capture_output=True, text=True)
    if result.returncode != 0:
        raise dodoException(
            f'PULCHRA failed with exit code {result.returncode}.\n'
            f'stdout:\n{result.stdout}\n'
            f'stderr:\n{result.stderr}'
        )

    if not rebuilt_path.is_file():
        raise dodoException(f'PULCHRA did not create expected output: {rebuilt_path}')

    if out_path is not None:
        Path(out_path).write_text(rebuilt_path.read_text())
        _normalize_protein_atom_elements(out_path)

    return rebuilt_path


def _ca_only_residue_indices(PDBParserObj):
    return [
        aa for aa in sorted(PDBParserObj.all_atom_coords_by_index)
        if set(PDBParserObj.all_atom_coords_by_index[aa]) == {'CA'}
    ]


def _merge_rebuilt_ca_only_residues(
    PDBParserObj,
    rebuilt_PDBParserObj,
    residue_indices,
    rebuilt_index_to_original_index,
):
    original_index_to_rebuilt_index = {
        original_index: rebuilt_index
        for rebuilt_index, original_index in rebuilt_index_to_original_index.items()
    }

    for aa in residue_indices:
        if aa not in original_index_to_rebuilt_index:
            raise dodoException(f'PULCHRA trace is missing residue index {aa}.')

        rebuilt_aa = original_index_to_rebuilt_index[aa]
        if rebuilt_aa not in rebuilt_PDBParserObj.all_atom_coords_by_index:
            raise dodoException(f'PULCHRA output is missing residue index {rebuilt_aa}.')
        if PDBParserObj.index_to_3aa[aa] != rebuilt_PDBParserObj.index_to_3aa[rebuilt_aa]:
            raise dodoException(
                f'Residue mismatch after PULCHRA at index {aa}: '
                f'{PDBParserObj.index_to_3aa[aa]} != {rebuilt_PDBParserObj.index_to_3aa[rebuilt_aa]}'
            )

        PDBParserObj.all_atom_coords_by_index[aa] = rebuilt_PDBParserObj.all_atom_coords_by_index[rebuilt_aa]
        PDBParserObj.all_atom_coords_by_index_with_aa[aa] = rebuilt_PDBParserObj.all_atom_coords_by_index_with_aa[rebuilt_aa]
        PDBParserObj.sequence_3aa_by_index[aa] = rebuilt_PDBParserObj.sequence_3aa_by_index[rebuilt_aa]

    PDBParserObj.number_atoms = sum(len(atoms) for atoms in PDBParserObj.all_atom_coords_by_index.values())
    return PDBParserObj


def write_ca_trace_from_PDBParserObj(PDBParserObj, out_path):
    residue_indices = sorted(PDBParserObj.all_atom_coords_by_index)
    ca_records = []

    for aa in residue_indices:
        atoms = PDBParserObj.all_atom_coords_by_index[aa]
        if 'CA' not in atoms:
            raise dodoException(f'Missing CA atom for residue index {aa}; cannot run PULCHRA.')
        ca_records.append({
            'xyz': atoms['CA'],
            'residue_name': PDBParserObj.index_to_3aa[aa],
            'beta': PDBParserObj.beta_vals_by_index[aa],
        })

    _write_ca_trace(ca_records, out_path)
    return {trace_index: residue_index for trace_index, residue_index in enumerate(residue_indices)}


def rebuild_PDBParserObj_with_pulchra(PDBParserObj, out_path, executable='pulchra',
    verbose=False, CONECT_lines=True):
    ca_only_residue_indices = _ca_only_residue_indices(PDBParserObj)
    if ca_only_residue_indices == []:
        raise dodoException('No CA-only residues were found for PULCHRA rebuilding.')

    with tempfile.TemporaryDirectory(prefix='dodo_pulchra_') as temp_dir:
        ca_trace_path = Path(temp_dir) / 'dodo_ca_trace.pdb'
        rebuilt_index_to_original_index = write_ca_trace_from_PDBParserObj(PDBParserObj, ca_trace_path)
        rebuilt_path = run_pulchra(ca_trace_path, executable=executable, verbose=verbose)
        rebuilt_PDBParserObj = PDBParser(rebuilt_path.read_text().splitlines())

    PDBParserObj = _merge_rebuilt_ca_only_residues(
        PDBParserObj,
        rebuilt_PDBParserObj,
        ca_only_residue_indices,
        rebuilt_index_to_original_index,
    )
    save_pdb_from_PDBParserObj(
        PDBParserObj,
        out_path=out_path,
        include_FD_atoms=True,
        CONECT_lines=CONECT_lines,
        add_mode='w',
        model_num=1,
        last_model=True,
    )
    _normalize_protein_atom_elements(out_path)


def rebuild_sequence_dict_with_pulchra(sequence_dict, out_path, executable='pulchra', verbose=False):
    with tempfile.TemporaryDirectory(prefix='dodo_pulchra_') as temp_dir:
        ca_trace_path = Path(temp_dir) / 'dodo_ca_trace.pdb'
        ca_records = [
            {'xyz': xyz, 'residue_name': residue_name, 'beta': 1.0}
            for xyz, residue_name in zip(sequence_dict['xyz_list'], sequence_dict['residue_names'])
        ]
        _write_ca_trace(ca_records, ca_trace_path)
        run_pulchra(ca_trace_path, out_path, executable=executable, verbose=verbose)
