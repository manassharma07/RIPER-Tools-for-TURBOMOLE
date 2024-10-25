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

# Set page config
st.set_page_config(page_title='`coord` File Visualizer', layout='wide', page_icon="⚛️",
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
            lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        
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
        return structure
    else:
        molecule = Molecule(atomic_species, coords)
        return molecule







# return filecontents
def read_file(filename):
    with open(filename, 'r') as file:
        return file.read()


# Streamlit app
st.write('# `coord` File Visualizer')
st.write(
    "#### Visualize Turbomole's `coord` files")  


st.write('You can either paste the `coord` file contents below or upload the `coord` file')
contents = st.text_area(label='Enter the contents of the `coord` file here', value='', placeholder='Put your text here',
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

    
    structure = parse_coord(contents)

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
        st.write('#### Fractional Coordinates of Center of Mass (Angstroms) ')
        st.write(com_fractional)

    