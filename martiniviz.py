#!/usr/bin/env python

"""
This script is designed for visualizing systems using the MARTINI v2 force field.
It reads a GROMACS *.top topology file to generate the bonds between beads and 
outputs a PSF file compatible with VMD for visualization.
"""

__author__ = "Ivo Kabelka"
__email__ = "ivo.kabelka@gmail.com"
__created__ = "2024-08-01"
__modified__ = "2024-08-01"
__version__ = "1.0.0"

import argparse
import sys
from collections import namedtuple

# Define a named tuple for storing atom information included in the PSF file
atom_tuple = namedtuple('atom', ('id', 'resid', 'resn', 'name', 'type'))

def read_top_file(top, included_files=set()):
    """
    Reads and processes a GROMACS topology file (*.top or *.itp).

    Args:
        top (file object): The topology file to read.
        included_files (set): Set of already included files to avoid circular dependencies.

    Returns:
        list: List of lines (without comments) from the topology file.
    """

    ifdefs = 0 # Counter for nested #ifdef blocks
    top_contents = [] # Concatenated contents of *.top and *.itp files

    for line in top:
        # Remove comments
        uncommented = line.split(';')[0].strip()

        # Skip empty lines
        if not uncommented:
            continue

        # Skip (nested) #ifdef blocks
        if uncommented.startswith(('#ifdef', '#ifndef')):
            ifdefs += 1
            continue
        elif uncommented.startswith('#endif'):
            ifdefs -= 1
            if ifdefs < 0:
                raise ValueError('Error: Invalid file format. Improper use of #ifdef / #endif blocks!')
            continue
        elif ifdefs > 0:
            continue

        # Include the contents of included *.itp files
        if uncommented.startswith('#include'):
            _, path = uncommented.split()
            path = path.strip('"\'')

            if path in included_files:
                raise ValueError(f'Error: Circular dependency detected: {path}')

            included_files.add(path)

            with open(path, 'r') as itp:
                # Recursively run this function to include *.itp file
                top_contents += read_top_file(itp)
        else:
            top_contents.append(uncommented)

    return top_contents

def split_section_molecules(uncommented):
    """
    Splits a line in the [ molecules ] section to extract the molecule name and count.

    Args:
        uncommented (str): The line to split.

    Returns:
        tuple: A tuple containing the molecule name and count.
    """

    try:
        mol_name, mol_count = uncommented.split()
        return mol_name, int(mol_count)
    except ValueError:
        print('Error: Invalid file format. Each line in the [ molecules ] section has to contain only the molecule_name and the molecule_count!')
        sys.exit(1)

def get_moleculetype_indices(topology):
    """
    Gets the indices of [ moleculetype ] definitions in the topology.

    Args:
        topology (list): The topology lines.

    Returns:
        list: Indices of moleculetype definitions.
    """

    indices = [e for e, row in enumerate(topology) if 'moleculetype' in row]
    indices.append(-1)

    return indices

def get_moleculetype_defines(topology):
    """
    Extracts bonded parameters from #define statements in the topology.

    Args:
        topology (list): The topology lines.

    Returns:
        dict: Dictionary of defines.
    """

    defines = {} # Dictionary of bonded parameters defined via a variable

    for row in topology:
        if row.startswith('#define'):
            r = row.split()
            if len(r) >= 3:
                defines[r[1]] = float(r[2])

    return defines

def parse_moleculetype_name(moltype):
    """
    Parses the name of the [ moleculetype ].

    Args:
        moltype (list): List of lines defining the moleculetype.

    Returns:
        str: Name of the moleculetype.
    """

    return moltype[1].split()[0]

def parse_section_name(row):
    """
    Parses the name of a [ section ].

    Args:
        row (str): The section line.

    Returns:
        str: Parsed section name.
    """

    return row[1:-1].strip().lower()

def split_atom_row(row):
    """
    Splits a row in the [ atoms ] section into its components.

    Args:
        row (str): The atom line.

    Returns:
        namedtuple: A namedtuple representing the atom.
    """

    a = row.split()

    return atom_tuple(int(a[0]), int(a[2]), a[3], a[4], a[1])

def split_bond_record(row, defines):
    """
    Splits a row in the [ bonds ] section into its components.

    Args:
        row (str): The bond line.
        defines (dict): Dictionary of defines for bond lengths.

    Returns:
        tuple: A tuple representing the bond.
    """

    b = row.split()

    # Bond or constrained length defined either in line or via a variable
    bond_length = defines[b[3]] if b[3] in defines.keys() else float(b[3])

    return int(b[0]), int(b[1]), bond_length

def parse_moleculetype(moltype, moltypes, defines):
    """
    Parses a [ moleculetype ] section.

    Args:
        moltype (list): List of lines defining the moleculetype.
        moltypes (dict): Dictionary to store parsed moleculetype data.
        defines (dict): Dictionary of defines for bond lengths.
    """

    # Get the name of the molecule
    mol_name = parse_moleculetype_name(moltype)

    moltypes[mol_name] = {'atoms': [], 'bonds': []} # Only atoms and bonds are needed in the PSF file
    bonds = [] # Temporary list of unsorted bonds
    section_name = '' # Name of the current section

    for row in moltype[2:]:
        if row[0] + row[-1] == '[]':
            section_name = parse_section_name(row)
        elif section_name == 'atoms':
            atom = split_atom_row(row)
            moltypes[mol_name]['atoms'].append(atom)
        elif section_name in ('bonds', 'constraints'):
            bond = split_bond_record(row, defines)
            if bond[2] < 0.5:
                bonds.append(bond[:2])

    # Sorted list of bonds; the dictionary is modified in-place and not returned
    moltypes[mol_name]['bonds'] = sorted(set(bonds), key=lambda x: (x[0], x[1]))

def parse_topology(topology):
    """
    Parses the entire topology to extract molecule types.

    Args:
        topology (list): List of topology lines.

    Returns:
        dict: Parsed molecule types.
    """

    # Get the indices of [ moleculetype ] sections in the topology list
    indices = get_moleculetype_indices(topology)

    # Get the bonded parameters defined via a variable
    defines = get_moleculetype_defines(topology)

    # Dictionary of all molecules
    moleculetypes = {}

    for i in range(1, len(indices)):
        moleculetype = topology[indices[i-1]:indices[i]]
        parse_moleculetype(moleculetype, moleculetypes, defines)

    return moleculetypes

def parse_section_molecules(molecules_raw):
    """
    Parses the [ molecules ] section to extract molecule names and counts.

    Args:
        molecules_raw (list): List of lines in the [ molecules ] section.

    Returns:
        list: List of tuples containing molecule names and counts.
    """

    molecules = []

    for mol in molecules_raw:
        molecules.append(split_section_molecules(mol))

    return molecules

def write_psf(out, molecules, moleculetypes):
    """
    Writes the PSF file with atoms and bonds information.

    Args:
        out (file object): Output file object.
        molecules (list): List of tuples containing molecule names and counts.
        moleculetypes (dict): Parsed molecule types.
    """

    atoms = [] # List of lines in the !NATOM header of the PSF file
    bonds = [] # LIST of lines following the !NBOND header in the PSF file

    atom_count = 1 # Atom counter used also for atom indices
    bond_count = 1 # Bond counter

    offset = 0 # Offset for bond indices to account for previous molecules

    # Iterate over all moleculetypes in the [ molecules ] section
    for mol_name, mol_count in molecules:
        # Append each moleculetype N times to reflect the number of coordinates in the system
        for _ in range(mol_count):
            for a in moleculetypes[mol_name]['atoms']:
                s = f'{atom_count:8d} {"P":4s} {a.resid:4d} {a.resn:4s} {a.name:4s} {a.type:41s}\n'
                atoms.append(s)
                atom_count += 1

            for b in moleculetypes[mol_name]['bonds']:
                s = f'{b[0] + offset:8d}{b[1] + offset:8d}'
                bonds.append(s)
                if bond_count % 4 == 0:
                    bonds.append('\n')
                bond_count += 1

            offset += len(moleculetypes[mol_name]['atoms'])

    # Write the PSF file
    out.write('PSF\n')
    out.write(f'{atom_count-1:8d} !NATOM\n')
    out.write(''.join(atoms))
    out.write(f'{bond_count-1:8d} !NBOND\n')
    out.write(''.join(bonds))

def parse_arguments():
    """
    Parses command-line arguments.

    Returns:
        Namespace: Parsed arguments.
    """

    parser = argparse.ArgumentParser(description='This script generates a PSF file with bonding information for the MARTINI force field.')
    parser.add_argument("-i", "--input", required=True, help="Specify the input TOP file.")
    parser.add_argument("-o", "--output", required=True, help="Specify the output PSF file.")

    return parser.parse_args()

def main():
    """
    Main function to read the input TOP file, process it, and write the output PSF file.
    """
    args = parse_arguments()

    with open(args.input, 'r') as top:
        top_contents = read_top_file(top)

    system_index = top_contents.index('[ system ]')
    topology = top_contents[:system_index]
    system = top_contents[system_index:]

    moleculetypes = parse_topology(topology)
    molecules_index = system.index('[ molecules ]')
    molecules = parse_section_molecules(system[molecules_index+1:])

    with open(args.output, 'w') as out:
        write_psf(out, molecules, moleculetypes)

if __name__ == '__main__':
    main()

