import streamlit as st
from mp_api.client import MPRester
from pymatgen.core import Structure, Element, Molecule, Lattice
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.cif import CifParser
from pymatgen.io.xyz import XYZ
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.pwscf import PWInput
import py3Dmol
import pandas as pd
import streamlit.components.v1 as components
from io import StringIO
from ase.io.espresso import read_espresso_in
from ase.io.extxyz import read_extxyz
from ase.io.cif import read_cif
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io.dmol import read_dmol_car
import re
import numpy as np
import io
# from mace.calculators import mace_mp

# Set page config
st.set_page_config(page_title='CIF/XYZ/CAR/POSCAR/PWSCF ➡️ RIPER', layout='wide', page_icon="⚛️",
                   menu_items={
                       'About': "A web app to help you with DFT related calculations using the RIPER module of [TURBOMOLE](https://www.turbomole.org/)"
                   })

# Sidebar stuff
st.sidebar.write('# About')
st.sidebar.write('Made By [Manas Sharma](https://manas.bragitoff.com)')
st.sidebar.write('In the group of [Prof. Dr. Marek Sierka](https://cmsg.uni-jena.de)')
st.sidebar.write('## Cite us:')
st.sidebar.write('[J. Phys. Chem. A 2025, 129, 39, 9062–9083](https://doi.org/10.1021/acs.jpca.5c02937)')
st.sidebar.write('### *Powered by*')
st.sidebar.write('* [Py3Dmol](https://3dmol.csb.pitt.edu/) for Chemical System Visualizations')
st.sidebar.write('* [Streamlit](https://streamlit.io/) for making of the Web App')
st.sidebar.write('* [PyMatgen](https://pymatgen.org/) for Periodic Structure Representations')
st.sidebar.write('* [PubChempy](https://pypi.org/project/PubChemPy/1.0/) for Accessing the PubChem Database')
st.sidebar.write('* [MP-API](https://pypi.org/project/mp-api/) for Accessing the Materials Project Database')
st.sidebar.write('* [ASE](https://wiki.fysik.dtu.dk/ase/) for File Format Conversions')
st.sidebar.write('### *Contributors*')
st.sidebar.write('[Ya-Fan Chen ](https://github.com/Lexachoc)')
st.sidebar.write('### *Source Code*')
st.sidebar.write('[GitHub Repository](https://github.com/manassharma07/RIPER-Tools-for-TURBOMOLE)')

@st.cache_resource
def get_mace_mp():
    return mace_mp(model="https://github.com/ACEsuit/mace-mp/releases/download/mace_omat_0/mace-omat-0-medium.model", device="cpu", default_dtype="float32")
    # return mace_mp(model="small", device="cpu", default_dtype="float32")

# Function to format floating-point numbers with alignment
def format_number(num, width=15, precision=8):
    # Handles positive/negative numbers while maintaining alignment
    # Adjusting the width based on the sign to ensure alignment
    return f"{num:>{width}.{precision}f}"


def calculate_com(structure):
    """
    Calculate the center of mass (COM) of a pymatgen structure.
    
    :param structure: pymatgen Structure object
    :return: COM coordinates as a tuple (x_com, y_com, z_com)
    """
    total_mass = 0
    com = [0, 0, 0]  # x_com, y_com, z_com
    
    for site in structure:
        element = Element(site.specie.symbol)
        mass = element.atomic_mass
        total_mass += mass
        
        # Multiply mass by the fractional coordinates
        com[0] += mass * site.coords[0]
        com[1] += mass * site.coords[1]
        com[2] += mass * site.coords[2]
    
    # Divide by total mass to get COM
    com = [x / total_mass for x in com]
    
    return tuple(com)

# Function to convert atomic coordinates to Bohr units
def convert_to_bohr(structure):
    coords = [(site.coords[0], site.coords[1], site.coords[2], site.species_string) for site in structure.sites]
    return [(x * 1.88972612456506, y * 1.88972612456506, z * 1.88972612456506, element.lower()) for x, y, z, element in coords]


# # Function to generate coordinate text
# def generate_coord_text(coords_bohr):
#     coord_text = "$coord\n"
#     for coord in coords_bohr:
#         coord_text += f"   {format_number(coord[0], precision=8)}  {format_number(coord[1], precision=8)}  {format_number(coord[2], precision=8)}  {coord[3]:<2s}\n"
#     coord_text += "$end"
#     return coord_text
# Function to generate coordinate text
def generate_coord_text(coords_bohr):
    coord_text = "$coord\n"
    for coord in coords_bohr:
        # Aligning all coordinates and ensuring the element symbol is aligned to the left
        coord_text += (
            f"{format_number(coord[0], width=15, precision=8)} "
            f"{format_number(coord[1], width=15, precision=8)} "
            f"{format_number(coord[2], width=15, precision=8)} "
            f"{coord[3]:<2s}\n"
        )
    coord_text += "$end"
    return coord_text


# Function to generate lattice parameter text
def generate_lattice_text(structure):
    lattice_params = structure.lattice.abc
    angles = structure.lattice.angles
    lattice_text = "$cell angs\n"
    lattice_text += f"  {lattice_params[0]:.8f}   {lattice_params[1]:.8f}   {lattice_params[2]:.8f}   {angles[0]}   {angles[1]}   {angles[2]}\n"
    lattice_text += "$periodic 3\n"
    lattice_text += "$kpoints\n"
    lattice_text += "    nkpoints 1 1 1 # Gamma point calculation"
    return lattice_text


# Function to convert a structure to CIF
def convert_to_cif(structure, filename):
    cif_writer = CifWriter(structure)
    cif_writer.write_file(filename)


# Function to visualize the structure using py3Dmol
def visualize_structure(structure, html_file_name='viz.html'):
    spin = st.checkbox('Spin', value=False, key='key' + html_file_name)
    view = py3Dmol.view(width=500, height=400)
    cif_for_visualization = structure.to(fmt="cif")
    view.addModel(cif_for_visualization, 'cif')
    # view.setStyle({'stick': {'radius': 0.2}})
    view.setStyle({'sphere': {'colorscheme': 'Jmol', 'scale': 0.3},
                   'stick': {'colorscheme': 'Jmol', 'radius': 0.2}})
    view.addUnitCell()
    view.zoomTo()
    view.spin(spin)
    view.setClickable({'clickable': 'true'})
    view.enableContextMenu({'contextMenuEnabled': 'true'})
    view.show()
    view.render()
    # view.png()
    t = view.js()
    f = open(html_file_name, 'w')
    f.write(t.startjs)
    f.write(t.endjs)
    f.close()

    HtmlFile = open(html_file_name, 'r', encoding='utf-8')
    source_code = HtmlFile.read()
    components.html(source_code, height=300, width=500)
    HtmlFile.close()


# Function to visualize the structure using py3Dmol
def visualize_molecule(structure, html_file_name='viz.html'):
    spin = st.checkbox('Spin', value=False, key='key' + html_file_name)
    view = py3Dmol.view(width=500, height=400)
    cif_for_visualization = structure.to(fmt="xyz")
    view.addModel(cif_for_visualization, 'xyz')
    # view.setStyle({'stick': {'radius': 0.2}})
    view.setStyle({'sphere': {'colorscheme': 'Jmol', 'scale': 0.3},
                   'stick': {'colorscheme': 'Jmol', 'radius': 0.2}})
    view.zoomTo()
    view.spin(spin)
    view.setClickable({'clickable': 'true'})
    view.enableContextMenu({'contextMenuEnabled': 'true'})
    view.show()
    view.render()
    # view.png()
    t = view.js()
    f = open(html_file_name, 'w')
    f.write(t.startjs)
    f.write(t.endjs)
    f.close()

    HtmlFile = open(html_file_name, 'r', encoding='utf-8')
    source_code = HtmlFile.read()
    components.html(source_code, height=300, width=500)
    HtmlFile.close()


# Function to display structure information
def display_structure_info(structure):
    st.subheader("Structure Information")
    st.write("Formula: ", structure.composition.reduced_formula)

    # Display lattice parameters
    a, b, c = structure.lattice.abc
    alpha, beta, gamma = structure.lattice.angles

    # Create a DataFrame for the lattice parameters and angles
    data = {
        "Lattice Parameters": [a, b, c, alpha, beta, gamma]
    }
    df_latt_params = pd.DataFrame(data, index=["a", "b", "c", "alpha", "beta", "gamma"])
    with st.expander("Lattice Parameters", expanded=False):
        st.table(df_latt_params)

    # Display lattice vectors
    lattice_vectors = structure.lattice.matrix
    df_vectors = pd.DataFrame(lattice_vectors, columns=["X", "Y", "Z"], index=["a", "b", "c"])
    with st.expander("Lattice Vectors", expanded=True):
        # st.write("Lattice Vectors:")
        st.table(df_vectors)

    # Create a list of atomic coordinates
    with st.expander("Atomic Coordinates", expanded=False):
        coord_type = st.selectbox('Coordinate type', ['Cartesian', 'Fractional/Crystal'])
        if coord_type == 'Cartesian':
            atomic_coords = []
            for site in structure.sites:
                atomic_coords.append([site.species_string] + list(site.coords))
        else:
            atomic_coords = []
            for site in structure.sites:
                atomic_coords.append([site.species_string] + list(site.frac_coords))

        # Create a Pandas DataFrame from the atomic coordinates list
        df_coords = pd.DataFrame(atomic_coords, columns=["Element", "X", "Y", "Z"])

        # Display the atomic coordinates as a table
        # st.write("Atomic Coordinates:")
        st.table(df_coords)


def build_lattice(a, b, c, alpha, beta, gamma):
    """Builds the lattice matrix from the cell parameters."""
    alpha, beta, gamma = np.radians([alpha, beta, gamma])
    volume = np.sqrt(1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 + 
                     2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma))
    lattice = np.zeros((3, 3))
    lattice[0] = [a, 0, 0]
    lattice[1] = [b * np.cos(gamma), b * np.sin(gamma), 0]
    lattice[2] = [c * np.cos(beta),
                  c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
                  c * volume / np.sin(gamma)]
    return lattice

def convert_pymatgen_to_ase_to_pymatgen(structure):
    convert_to_cif(structure, "temp.cif")
    file = open("temp.cif", 'r')
    return parse_cif_ase(file)

def parse_coord(file_contents):
    lines = file_contents.strip().splitlines()

    lattice = None
    coords = []
    atomic_species = []
    fractional = False
    in_bohr = False
    is_periodic = False

    for i, line in enumerate(lines):
        line = line.strip()
        
        # Check for units and periodicity
        if line.startswith('$cell'):
            is_periodic = True
            if 'angs' in line:
                in_bohr = False
            else:
                in_bohr = True
            
            # Extract cell parameters (axes and angles)
            cell_line = lines[i + 1].strip().split()
            a, b, c = map(float, cell_line[:3])
            alpha, beta, gamma = map(float, cell_line[3:])
            # lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
            lattice = build_lattice(a, b, c, alpha, beta, gamma)
        
        elif line.startswith('$lattice'):
            is_periodic = True
            if 'angs' in line:
                in_bohr = False
            else:
                in_bohr = True

            # Parse lattice vectors from following 3 lines
            vec1 = list(map(float, lines[i + 1].strip().split()))
            vec2 = list(map(float, lines[i + 2].strip().split()))
            vec3 = list(map(float, lines[i + 3].strip().split()))
            lattice = Lattice([vec1, vec2, vec3])
        
        elif line.startswith('$coord frac'):
            is_periodic = True
            fractional = True
        
        elif line.startswith('$coord'):
            # Parse coordinates (in Cartesian or fractional)
            for j in range(i + 1, len(lines)):
                coord_line = lines[j].strip()
                if coord_line == '$end':
                    break
                parts = coord_line.split()
                
                if len(parts) == 4:
                    x, y, z = map(float, parts[:3])
                    atomic_species.append(parts[3].capitalize())
                    coords.append([x, y, z])
        
    # Convert Bohr to Angstrom if needed
    if not fractional:
        coords = [[x * 0.52917721092, y * 0.52917721092, z * 0.52917721092] for x, y, z in coords]

    if in_bohr:
        lattice = Lattice(lattice.matrix * 0.52917721092)
    
    # Create Structure or Molecule
    if is_periodic:
        structure = Structure(lattice, atomic_species, coords, coords_are_cartesian=not fractional)
        structure = convert_pymatgen_to_ase_to_pymatgen(structure)
        return structure
    else:
        molecule = Molecule(atomic_species, coords)
        return molecule

def parse_cif_pymatgen(contents):
    # Parse the CIF file using pymatgen
    cif_parser = CifParser.from_str(contents)
    structure = cif_parser.get_structures(primitive=False)[0]  # Assuming there's only one structure in the CIF file
    return structure


def parse_xyz(contents):
    # Parse the XYZ file using pymatgen
    xyz_parser = XYZ.from_str(contents)
    structure = xyz_parser.molecule  # Assuming it's a molecule XYZ file
    return structure


def parse_poscar(contents):
    # Parse the POSCAR file using pymatgen
    poscar_parser = Poscar.from_str(contents)
    structure = poscar_parser.structure
    return structure


def parse_quantum_espresso(contents):
    # Parse the POSCAR file using pymatgen
    qe_parser = PWInput.from_str(contents)
    structure = qe_parser.structure
    return structure


def parse_qe_ase(stringio):
    # Read QE input file
    atoms = read_espresso_in(stringio)

    # Convert ASE Atoms to pymatgen Structure
    structure = AseAtomsAdaptor().get_structure(atoms)
    return structure

def parse_extxyz_ase(stringio):
    # Read extended XYZ file
    # atoms = read_extxyz(stringio)
    atoms = read(stringio, format="extxyz")

    # Convert ASE Atoms to pymatgen Structure
    structure = AseAtomsAdaptor().get_structure(atoms)
    return structure


def parse_cif_ase(stringio):
    # Read CIF
    atoms = read(stringio, format="cif")

    # Convert ASE Atoms to pymatgen Structure
    structure = AseAtomsAdaptor().get_structure(atoms)
    return structure

def parse_car_ase(stringio):
    """Parse a CAR file, handling PBC=2D, PBC=3D, and non-periodic cases."""
    
    # Read the content to detect PBC type
    content = stringio.read()
    stringio.seek(0)  # Reset for reading
    
    # Detect PBC dimensionality
    pbc_match = re.search(r'PBC=(\S+)', content)
    pbc_type = pbc_match.group(1) if pbc_match else None
    
    if pbc_type == '2D':
        # For 2D PBC, we need to modify the file to look like 3D so ASE can parse it
        # Replace PBC=2D with PBC=ON so ASE's parser can handle it
        modified_content = content.replace('PBC=2D', 'PBC=ON')
        modified_stringio = io.StringIO(modified_content)
        
        atoms = read_dmol_car(modified_stringio)
        
        # Parse the PBC line to get 2D lattice parameters
        # Format: PBC    a    b    angle (space_group)
        pbc_line = None
        for line in content.splitlines():
            if line.startswith('PBC') and not line.startswith('PBC='):
                pbc_line = line
                break
        
        if pbc_line is not None:
            # Extract lattice parameters from PBC line
            # For 2D: PBC  a  b  gamma(space_group)
            parts = pbc_line.split()
            # parts[0] = 'PBC', then a, b, gamma(sg)
            a = float(parts[1])
            b = float(parts[2])
            # The angle might have the space group appended like "93.0000(p 1)"
            angle_str = parts[3]
            gamma = float(re.match(r'([0-9.]+)', angle_str).group(1))
            
            # Set 2D periodicity (periodic in a, b; non-periodic in c)
            import numpy as np
            
            # Build 2D lattice vectors in xy-plane
            ax = a
            ay = 0.0
            bx = b * np.cos(np.radians(gamma))
            by = b * np.sin(np.radians(gamma))
            
            # Set a large vacuum in z direction
            # Find the range of z coordinates to set appropriate vacuum
            z_coords = atoms.positions[:, 2]
            z_range = z_coords.max() - z_coords.min()
            c_vacuum = max(20.0, z_range + 15.0)  # At least 20 Å or z_range + 15 Å
            
            cell = np.array([
                [ax, ay, 0.0],
                [bx, by, 0.0],
                [0.0, 0.0, c_vacuum]
            ])
            
            atoms.set_cell(cell)
            atoms.set_pbc([True, True, False])
        
        # Convert to pymatgen Structure with 2D periodicity
        # Since pymatgen Structure requires 3D periodicity, we treat it as a slab
        structure = AseAtomsAdaptor().get_structure(atoms)
        return structure
    
    else:
        # Standard 3D or non-periodic case
        atoms = read_dmol_car(stringio)
        
        if any(atoms.pbc):
            structure = AseAtomsAdaptor().get_structure(atoms)
        else:
            structure = AseAtomsAdaptor().get_molecule(atoms)
        return structure

# def parse_car_ase(stringio):
#     # Read CAR
#     atoms = read_dmol_car_custom(stringio)

#     # Convert ASE Atoms to pymatgen Structure (determine if CAR file is 3D periodicity or not)
#     if any(atoms.pbc):
#         structure = AseAtomsAdaptor().get_structure(atoms)
#     else:
#         structure = AseAtomsAdaptor().get_molecule(atoms)
#     return structure

def read_dmol_car_custom(stringio):
    """
    Parse a .car file and return an ASE Atoms object supporting 0, 1, 2, or 3D PBC.
    
    ASE's read_dmol_car only supports 000 or 111 PBC. This custom parser
    reads the PBC and cell vectors directly from the file header.
    """
    import numpy as np
    from ase import Atoms

    lines = stringio.read().splitlines()

    # --- Parse header ---
    pbc = [False, False, False]
    cell = None
    periodicity_map = {
        "cluster":    [False, False, False],
        "polymer":    [True,  False, False],
        "slab":       [True,  True,  False],
        "crystal":    [True,  True,  True],
    }

    i = 0
    for i, line in enumerate(lines):
        stripped = line.strip().lower()

        # Periodicity keyword (line starting with "!DATE" marks end of header)
        for keyword, pbc_val in periodicity_map.items():
            if stripped.startswith(keyword):
                pbc = pbc_val
                break

        # Cell vectors: lines starting with "PBC" contain a, b, c, alpha, beta, gamma
        if stripped.startswith("pbc") and not stripped.startswith("pbc="):
            parts = line.split()
            if len(parts) >= 7:
                try:
                    a, b, c = float(parts[1]), float(parts[2]), float(parts[3])
                    alpha, beta, gamma = float(parts[4]), float(parts[5]), float(parts[6])
                    from ase.geometry import cellpar_to_cell
                    cell = cellpar_to_cell([a, b, c, alpha, beta, gamma])
                except ValueError:
                    pass

        if stripped.startswith("!date"):
            i += 1  # atom lines start after this
            break

    # --- Parse atom lines ---
    symbols = []
    positions = []

    for line in lines[i:]:
        stripped = line.strip()
        if not stripped or stripped.lower().startswith("end"):
            break
        parts = stripped.split()
        # CAR atom line format: name x y z mol_name resname resnum element charge
        if len(parts) < 8:
            continue
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            continue

        # Element symbol is in column index 7
        element = parts[7]
        # Strip trailing digits (e.g. "C1" -> "C")
        element = ''.join(c for c in element if c.isalpha())
        symbols.append(element)
        positions.append([x, y, z])

    atoms = Atoms(
        symbols=symbols,
        positions=positions,
        cell=cell if cell is not None else np.zeros((3, 3)),
        pbc=pbc,
    )

    return atoms


# return filecontents
def read_file(filename):
    with open(filename, 'r') as file:
        return file.read()


# Streamlit app
st.write('# CIF/XYZ/CAR/POSCAR/PWSCF ➡️ RIPER')
st.write(
    "#### Get atomic coordinates and cell parameters for RIPER (TURBOMOLE) from a CIF/POSCAR or XYZ (only atomic coordinates)")  # TODO

st.write("Please select the file format")

# Select file format
file_format = st.selectbox("Select file format",
                           ("CIF", "XYZ", "CAR (Materials Studio)", "POSCAR", "Quantum ESPRESSO (PWSCF)", "Extended XYZ", "TMOL (Coord)"))

if file_format =="TMOL (Coord)":
    st.warning('Only molecular or coord files with 3D periodicity are supported. That is, the cell parameters of a 3D cell can only be parsed.')

if file_format == 'CIF':
    cif_parser_options = ['ASE', 'PYMATGEN']
    selected_cif_parser = st.selectbox("Select a parser for CIFs", cif_parser_options)

st.write('You can either paste the source file contents below or upload the source file')
contents = st.text_area(label='Enter the contents of the source file here', value='', placeholder='Put your text here',
                        height=400, key='input_text_area')

# Create a file uploader widget
file = st.file_uploader("or Upload the file")

if file is not None:
    # If a file is uploaded, read its contents
    # contents = file.read()
    # To read file as bytes:
    bytes_data = file.getvalue()

    # To convert to a string based IO:
    stringio = StringIO(file.getvalue().decode("utf-8"))

    # To read file as string:
    contents = stringio.read()

if contents != '':

    # Parse the file based on the selected format
    if file_format == "CIF":
        if selected_cif_parser == 'PYMATGEN':
            structure = parse_cif_pymatgen(contents)
        elif selected_cif_parser == 'ASE':
            # Create a StringIO object
            stringio_obj_cif = StringIO(contents)
            structure = parse_cif_ase(stringio_obj_cif)
    elif file_format == "XYZ":
        structure = parse_xyz(contents)
    elif file_format == "CAR (Materials Studio)":
        stringio_obj_car = StringIO(contents)
        # try:
        structure = parse_car_ase(stringio_obj_car)
        # except Exception:
        #     raise Exception("Wrong CAR format or PBC is not set to ON or OFF (PBC=2D is not supported here).")
    elif file_format == "POSCAR":
        structure = parse_poscar(contents)
    elif file_format == "Quantum ESPRESSO (PWSCF)":
        # structure = parse_quantum_espresso(contents)
        # Create a StringIO object
        stringio_obj = StringIO(contents)
        structure = parse_qe_ase(stringio_obj)
    elif file_format == "Extended XYZ":
        # Create a StringIO object
        stringio_obj = StringIO(contents)
        structure = parse_extxyz_ase(stringio_obj)
    elif file_format == "TMOL (Coord)":
        structure = parse_coord(contents)

    # if file_format!="XYZ" and selected_cif_parser=='PYMATGEN':
    #     # Get conventional structure
    #     analyzer = SpacegroupAnalyzer(structure)
    #     structure = analyzer.get_conventional_standard_structure()

    st.write("Successfully parsed file!")
    # st.write("Pymatgen Structure:")
    # st.write(structure)

    # Display structure information
    if isinstance(structure, Structure):  # type = Structure
        display_structure_info(structure)
    else:  # type = Molecule
        st.subheader("Atomic Coordinates")
    # if not file_format == 'XYZ' :
    #     display_structure_info(structure)
    # elif file_format == 'CAR (Materials Studio)':
    #     if isinstance(structure, Structure):  # type = Structure
    #         display_structure_info(structure)
    #     else:  # type = Molecule
    #         st.subheader("Atomic Coordinates")
    # else:
    #     st.subheader("Atomic Coordinates")
    #     # Create a dataframe with atomic symbols and atomic coordinates

    # Visualize the structure
    if isinstance(structure, Structure):
        visualize_structure(structure, "viz1.html")
    else:  # type = Molecule
        visualize_molecule(structure, "viz1.html")

    # Download CIF files
    if isinstance(structure, Structure):
        st.subheader("Download CIF File")

        convert_to_cif(structure, "structure.cif")
        st.download_button('Download CIF', data=read_file("structure.cif"), file_name='structure.cif', key='cif_button')

    center_of_mass = calculate_com(structure)
    st.write('#### Cartesian Coordinates of Center of Mass (Angstroms) ')
    st.write(center_of_mass)
    if isinstance(structure, Structure):
        # Convert COM in Cartesian coordinates to fractional coordinates
        com_fractional = structure.lattice.get_fractional_coords(center_of_mass)
        st.write('#### Fractional Coordinates of Center of Mass ')
        st.write(com_fractional)

    # Get TURBOMOLE (RIPER) Coord file and Control file contents
    st.subheader("RIPER Files")
    # Convert the atomic coordinates to Bohr units
    coords_bohr = convert_to_bohr(structure)

    # Generate the coordinate text
    coords_text = generate_coord_text(coords_bohr)

    # Generate the lattice parameter text
    if isinstance(structure, Structure):
        lattice_text = generate_lattice_text(structure)

    # Create two columns for text boxes
    col1, col2 = st.columns(2)

    # Display the coordinate text in the first column
    with col1:
        st.text_area("Coord file contents (Cartesian coordinates in Bohr)", value=coords_text, height=300,
                     key='coords_text')
        st.download_button('Download coord file', coords_text, file_name='coord', key='control_text')

    if isinstance(structure, Structure):
        # Display the lattice parameters text in the second column
        with col2:
            st.text_area("Add the following to your control file", value=lattice_text, height=300)
    
    # if st.button("Calculate energy with MACE (ML) model trained on OMAT Dataset"):
    #     mace_mp = get_mace_mp()
    #     # Convert pymatgen structure to ASE Atoms object and evaluate energy using mace model
    #     structure_ase_temp = AseAtomsAdaptor().get_atoms(structure)
    #     structure_ase_temp.set_calculator(mace_mp)
    #     energy = structure_ase_temp.get_potential_energy()
    #     st.write("Energy (eV): ", energy)
    #     # write forces as a df table
    #     forces = structure_ase_temp.get_forces()
    #     st.write("Forces (eV/Ang): ")
    #     st.write(forces)

    if isinstance(structure, Structure):
        # Add sliders for translation along a, b, and c lattice vectors
        # st.warning('Translation feature is still in development so may not work as expected!')
        st.subheader("Translate Structure Along Lattice Vectors")
        translate_a = st.slider("Translate along a", min_value=-1.0, max_value=1.0, step=0.01, value=0.0)
        translate_b = st.slider("Translate along b", min_value=-1.0, max_value=1.0, step=0.01, value=0.0)
        translate_c = st.slider("Translate along c", min_value=-1.0, max_value=1.0, step=0.01, value=0.0)

        translated_structure = structure.copy()
        # Apply the translation to the structure
        translation_vector = translate_a * structure.lattice.matrix[0] + \
                            translate_b * structure.lattice.matrix[1] + \
                            translate_c * structure.lattice.matrix[2]
        translated_structure.translate_sites(range(len(structure.sites)), translation_vector, frac_coords=False)

        # Re-visualize the translated structure

        # Get the number of atoms
        num_atoms_supercell = translated_structure.num_sites
        if num_atoms_supercell<500:
            visualize_structure(translated_structure, 'viz_translated.html')
        else:
            st.warning("We can't visualize your cell as it contains more than 500 atoms which is a bit too much for a free web app.\n But don't worry, RIPER can still do the calculations with ease (provided you have the required resources).")

        # Download CIF files
        if isinstance(translated_structure, Structure):
            st.subheader("Download CIF File for the translated structure")

            convert_to_cif(translated_structure, "translated_structure.cif")
            st.download_button('Download CIF', data=read_file("translated_structure.cif"), file_name='translated_structure.cif', key='cif_button_translated_structure')

        # Get TURBOMOLE (RIPER) Coord file and Control file contents
        st.subheader("RIPER Files for the Translated Structure")
        # Convert the atomic coordinates to Bohr units
        coords_bohr_translated = convert_to_bohr(translated_structure)

        # Generate the coordinate text
        coords_text_translated = generate_coord_text(coords_bohr_translated)
        # Generate the lattice parameter text
        lattice_text_translated = generate_lattice_text(translated_structure)

        # Create two columns for text boxes
        col1, col2 = st.columns(2)

        # Display the coordinate text in the first column
        with col1:
            st.text_area("Coord file contents (Translated Cartesian coordinates in Bohr)", value=coords_text_translated, height=300)

            st.download_button('Download translated coord file', coords_text_translated, file_name='coord')

        # Display the lattice parameters text in the second column
        with col2:
            st.text_area("Add the following to your control file (Translated)", value=lattice_text_translated, height=300)
        


    if isinstance(structure, Structure):
        # Create supercells
        # with st.expander("Model Supercell", expanded=False):
        st.subheader('Model Supercell')
        # Create three columns for inputs
        col1, col2, col3 = st.columns(3)

        # Add input fields for nx, ny, and nz in separate columns
        nx = col1.number_input("Enter nx", min_value=1, step=1)
        ny = col2.number_input("Enter ny", min_value=1, step=1)
        nz = col3.number_input("Enter nz", min_value=1, step=1)
        # supercell_structure = structure.copy()
        supercell_structure = translated_structure.copy()
        supercell_structure.make_supercell([int(nx), int(ny), int(nz)])
        # Get the number of atoms
        num_atoms_supercell = supercell_structure.num_sites
        if num_atoms_supercell<500:
            visualize_structure(supercell_structure, 'viz2.html')
        else:
            st.warning("We can't visualize your supercell as it contains more than 500 atoms which is a bit too much for a free web app.\n But don't worry, RIPER can still do the calculations with ease (provided you have the required resources).")


        # Get TURBOMOLE (RIPER) Coord file and Control file contents
        st.subheader("RIPER Files for the Supercell")
        # Convert the atomic coordinates to Bohr units
        coords_bohr_super = convert_to_bohr(supercell_structure)

        # Generate the coordinate text
        coords_text_super = generate_coord_text(coords_bohr_super)
        # Generate the lattice parameter text
        lattice_text_super = generate_lattice_text(supercell_structure)

        # Create two columns for text boxes
        col1, col2 = st.columns(2)

        # Display the coordinate text in the first column
        with col1:
            st.text_area("Coord file contents (Supercell Cartesian coordinates in Bohr)", value=coords_text_super, height=300)

            st.download_button('Download supercell coord file', coords_text_super, file_name='coord')

        # Display the lattice parameters text in the second column
        with col2:
            st.text_area("Add the following to your control file (Supercell)", value=lattice_text_super, height=300)
