from ase import Atoms
from ase.build import make_supercell
from abtem.structures import orthogonalize_cell, is_cell_orthogonal
from ase.io import read, write
import numpy as np

def remove_overlapping_atoms(atoms, tolerance=0.1):
    """Remove overlapping atoms within a certain tolerance."""
    positions = atoms.get_positions()
    unique_indices = []
    for i, pos1 in enumerate(positions):
        if all(np.linalg.norm(pos1 - positions[j]) > tolerance for j in unique_indices):
            unique_indices.append(i)
    return atoms[unique_indices]

def main(input_cif, output_cif, supercell_matrix=None):
    strain = None
    # Read the initial structure from a CIF file
    atoms = read(input_cif)

    # Apply the supercell transformation if a supercell matrix is provided
    if supercell_matrix is not None:
        atoms = make_supercell(atoms, supercell_matrix)
        print(f"Supercell created with matrix:\n{supercell_matrix}")

    # Check if the cell is already orthogonal
    if is_cell_orthogonal(atoms):
        print("The cell is already orthogonal.")
        orthogonal_atoms = atoms
        strain = None
    else:
        # Perform the conversion to an orthogonal cell
        orthogonal_atoms, strain = orthogonalize_cell(atoms, max_repetitions=5, return_transform=True)
        # orthogonal_atoms = orthogonalize_cell(atoms)


    # Print the strain and new cell parameters if applicable
    if strain is not None:
        print('Strain:\n', strain)
    print('Orthogonalized Cell:\n', orthogonal_atoms)
    print("Transformed coordinates:")
    print(orthogonal_atoms.get_positions())
    print("Transformed cell:")
    print(orthogonal_atoms.get_cell())

    # Save the orthogonalized structure to a CIF file
    write(output_cif, orthogonal_atoms)

if __name__ == "__main__":
    # Example usage
    input_cif = 'input.cif'  # Replace with your CIF file name
    output_cif = 'orthogonalized_structure.cif'  # Replace with desired output file name

    # Define a supercell matrix if needed (e.g., 2x2x1 supercell)
    supercell_matrix = [[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, 1]]

    main(input_cif, output_cif, supercell_matrix)
