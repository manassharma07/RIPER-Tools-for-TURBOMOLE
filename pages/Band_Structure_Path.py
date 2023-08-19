import streamlit as st
import numpy as np
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor


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

def create_structure_from_parameters(a, b, c, alpha, beta, gamma):
    a_vec = [a, 0, 0]
    b_vec = [b * np.cos(np.radians(gamma)), b * np.sin(np.radians(gamma)), 0]
    c_x = c * np.cos(np.radians(beta))
    c_y = c * (np.cos(np.radians(alpha)) - np.cos(np.radians(beta)) * np.cos(np.radians(gamma))) / np.sin(np.radians(gamma))
    c_z = np.sqrt(c**2 - c_x**2 - c_y**2)
    c_vec = [c_x, c_y, c_z]
    
    lattice_matrix = np.array([a_vec, b_vec, c_vec])
    lattice = lattice_matrix.T
    
    return Structure(lattice, ["H"], [[0, 0, 0]])

# def generate_turbomole_text(bandpath_str):
#     special_points = bandpath_str.split(",")
#     kptlines = len(special_points) - 1
#     text = f"kptlines {kptlines}\n"
    
#     for i in range(len(special_points) - 1):
#         start_point = special_points[i]
#         end_point = special_points[i + 1]
#         line = f"recipr {start_point} {end_point} 40\n"
#         text += line
    
#     return text
def calculate_nlines(bandpath_str):
    substrings = bandpath_str.split(",")
    nlines = sum(len(substring) - 1 for substring in substrings)
    return nlines

def generate_turbomole_text(bandpath_str, special_points):
    substrings = bandpath_str.split(",")
    text = ""

    for substring in substrings:
        n_chars = len(substring)

        for i in range(n_chars - 1):
            start_point = special_points[substring[i]]
            end_point = special_points[substring[i + 1]]
            line = f"recipr {start_point} {end_point} 40\n"
            text += line

    text = f"kptlines {calculate_nlines(bandpath_str)}\n" + text
    return text

st.title('Band Structure Path')

structure = None

cif_file = st.file_uploader("Upload CIF file")
cif_contents = st.text_area("Or paste CIF contents here")

use_primitive = st.checkbox("Convert to primitive cell?", value=True)
if not use_primitive:
    st.warning("Warning: Band structure might be affected by band folding if not using primitive cell.")

lattice_parameters_manual = st.checkbox("Enter lattice parameters manually")
if lattice_parameters_manual:
    a = st.number_input("a", value=1.0)
    b = st.number_input("b", value=1.0)
    c = st.number_input("c", value=1.0)
    alpha = st.number_input("alpha", value=90.0)
    beta = st.number_input("beta", value=90.0)
    gamma = st.number_input("gamma", value=90.0)
    structure = create_structure_from_parameters(a, b, c, alpha, beta, gamma)
else:
    if cif_file:
        cif_text = cif_file.read().decode("utf-8")
        structure = Structure.from_str(cif_text, fmt="cif")
    elif cif_contents:
        structure = Structure.from_str(cif_contents, fmt="cif")
    else:
        st.write("Please upload a CIF file or paste its contents.")

if structure:
    if use_primitive:
        primitive_structure = structure.get_primitive_structure()
        st.success("Converted to Primitive Structure! Using primitive structure from now on.")
        atoms = AseAtomsAdaptor.get_atoms(primitive_structure)
    else:
        st.warning("Using Conventional Structure. May result in Band Folding")
        atoms = AseAtomsAdaptor.get_atoms(structure)

    lat = atoms.cell.get_bravais_lattice()
    special_points = lat.get_special_points()

    st.write("### Special High-Symmetry k-points:")
    st.write(special_points)

    bandpath = atoms.cell.bandpath()
    bandpath_str = bandpath.path
    st.write("### Special High-Symmetry path for Band Structure Calculation:")
    st.write(bandpath_str)
    nlines = calculate_nlines(bandpath_str)

    st.write("#### Number of lines (paths) in band structure:", nlines)

    st.write("### Input text for RIPER band structure calculation (Add it to your `control` file)")
    text_area_content = generate_turbomole_text(bandpath_str)
    turbomole_text = st.text_area("TURBOMOLE input text", value=text_area_content, height=200)



