import streamlit as st
from mp_api.client import MPRester
from pymatgen.core import Structure
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

# Set page config
st.set_page_config(page_title='CIF/XYZ/CAR/POSCAR/PWSCF ➡️ RIPER', layout='wide', page_icon="⚛️",
                   menu_items={
                       'About': "A web app to help you with DFT related calculations using the RIPER module of [TURBOMOLE](https://www.turbomole.org/)"
                   })

# Sidebar stuff
st.sidebar.write('# About')
st.sidebar.write('Made By [Manas Sharma](https://manas.bragitoff.com)')
st.sidebar.write('In the group of [Prof. Dr. Marek Sierka](https://cmsg.uni-jena.de)')
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


# Function to convert atomic coordinates to Bohr units
def convert_to_bohr(structure):
    coords = [(site.coords[0], site.coords[1], site.coords[2], site.species_string) for site in structure.sites]
    return [(x * 1.88972612456506, y * 1.88972612456506, z * 1.88972612456506, element.lower()) for x, y, z, element in coords]


# Function to generate coordinate text
def generate_coord_text(coords_bohr):
    coord_text = "$coord\n"
    for coord in coords_bohr:
        coord_text += f"{coord[0]:>20.14f}  {coord[1]:>20.14f}  {coord[2]:>20.14f}  {coord[3]:<2s}\n"
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
    components.html(source_code, height=300, width=900)
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
    components.html(source_code, height=300, width=900)
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
    atoms = read(stringio, format="read_extxyz")

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
    # Read CAR
    atoms = read_dmol_car(stringio)

    # Convert ASE Atoms to pymatgen Structure (determine if CAR file is 3D periodicity or not)
    if all(atoms.pbc):
        structure = AseAtomsAdaptor().get_structure(atoms)
    else:
        structure = AseAtomsAdaptor().get_molecule(atoms)
    return structure



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
                           ("CIF", "XYZ", "CAR (Materials Studio)", "POSCAR", "Quantum ESPRESSO (PWSCF)", "Extended XYZ"))

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
        try:
            structure = parse_car_ase(stringio_obj_car)
        except Exception:
            raise Exception("Wrong CAR format or PBC is not set to ON or OFF (PBC=2D is not supported here).")
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
        structure = parse_extxyz_ase(contents)

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
        st.subheader("Download CIF Files")

        convert_to_cif(structure, "structure.cif")
        st.download_button('Download CIF', data=read_file("structure.cif"), file_name='structure.cif', key='cif_button')

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
        # st.write(primitive_structure)
        # st.write(conventional_structure)
        supercell_structure = structure.copy()
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
            st.text_area("Coord file contents (Cartesian coordinates in Bohr)", value=coords_text_super, height=300,
                         key='supercell_text_coord')

            st.download_button('Download coord file', coords_text_super, file_name='coord')

        # Display the lattice parameters text in the second column
        with col2:
            st.text_area("Add the following to your control file", value=lattice_text_super, height=300,
                         key='supercell_text_control')
