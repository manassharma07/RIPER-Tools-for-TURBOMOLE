from ase.lattice.cubic import SimpleCubic, FaceCenteredCubic, BodyCenteredCubic
from ase.lattice.tetragonal import SimpleTetragonal, CenteredTetragonal
from ase.lattice.orthorhombic import (SimpleOrthorhombic, BaseCenteredOrthorhombic, 
                                      FaceCenteredOrthorhombic, BodyCenteredOrthorhombic)
from ase.lattice.monoclinic import SimpleMonoclinic, BaseCenteredMonoclinic
from ase.lattice.triclinic import Triclinic
from ase.lattice.hexagonal import Hexagonal
# from ase.dft.kpoints import bandpath
from ase.lattice import BCT, CUB

# def print_high_symmetry_points_and_bandpath(lattice_name, lattice):
#     print(BCT(3, 5).bandpath())
#     print(lattice.lattice)
#     special_points = lattice.get_special_points()
    
#     # Get the band path
#     bandpath_obj = lattice.bandpath()
#     bandpath_str = lattice.bandpath().path

#     print(f"\nBravais Lattice Type: {lattice_name}")
#     print(f"High-Symmetry k-points: {special_points}")
#     print(f"Band Path: {bandpath_str}")
#     print("="*50)
def print_high_symmetry_points_and_bandpath(lattice_name, lattice):
    # Generate the ASE atoms object from the lattice
    atoms = lattice#.atoms()

    # Get the Bravais lattice and high-symmetry k-points
    bravais_lattice = atoms.cell.get_bravais_lattice(pbc=[1, 1, 1])
    special_points = bravais_lattice.get_special_points()
    
    # Get the band path
    bandpath_obj = atoms.cell.bandpath()
    bandpath_str = bandpath_obj.path

    print(f"\nBravais Lattice Type: {lattice_name}")
    print(f"High-Symmetry k-points: {special_points}")
    print(f"Band Path: {bandpath_str}")
    print("="*50)

def generate_bravais_lattice_data():
    # Define the 14 3D Bravais lattices using the ase.lattice module
    bravais_lattices = {
        "Simple Cubic (SC)": SimpleCubic(latticeconstant=3.6, symbol='Cu'),
        "Face-Centered Cubic (FCC)": FaceCenteredCubic(latticeconstant=3.6, symbol='Cu'),
        "Body-Centered Cubic (BCC)": BodyCenteredCubic(latticeconstant=2.87, symbol='Cu'),
        "Simple Tetragonal": SimpleTetragonal(latticeconstant=[3.5, 5.6], symbol='Cu'),
        "Body-Centered Tetragonal (BCT)": CenteredTetragonal(latticeconstant=[3.5, 5.6], symbol='Cu'),
        # "Simple Orthorhombic": SimpleOrthorhombic(a=5.2, b=6.3, c=7.1, symbol='Cu'),
        # "Base-Centered Orthorhombic": BaseCenteredOrthorhombic(a=4.2, b=5.6, c=7.9, symbol='Cu'),
        # "Face-Centered Orthorhombic": FaceCenteredOrthorhombic(a=4.9, b=6.7, c=8.1, symbol='Cu'),
        # "Body-Centered Orthorhombic": BodyCenteredOrthorhombic(a=3.9, b=4.9, c=6.7, symbol='Cu'),
        # "Simple Monoclinic": SimpleMonoclinic(a=4.5, b=5.2, c=5.4, alpha=90, beta=110, gamma=90, symbol='Cu'),
        # "Base-Centered Monoclinic": BaseCenteredMonoclinic(a=4.2, b=6.2, c=7.3, alpha=90, beta=103, gamma=90, symbol='Cu'),
        # "Triclinic": Triclinic(a=4.3, b=5.1, c=6.2, alpha=60, beta=70, gamma=80, symbol='Cu'),
        # "Hexagonal": Hexagonal(a=3.21, c=5.21, symbol='Cu'),
    }

    # Loop through each Bravais lattice, printing high-symmetry points and band paths
    for lattice_name, lattice in bravais_lattices.items():
        print_high_symmetry_points_and_bandpath(lattice_name, lattice)

if __name__ == "__main__":
    generate_bravais_lattice_data()
