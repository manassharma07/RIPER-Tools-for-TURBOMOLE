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


def calculate_nlines(bandpath_str):
    substrings = bandpath_str.split(",")
    nlines = sum(len(substring) - 1 for substring in substrings)
    return nlines

def format_coordinate(coord):
    return "    ".join(f"{c:.6f}" for c in coord)

def generate_turbomole_text(nlines, bandpath_str, special_points):
    substrings = bandpath_str.split(",")
    text = ""

    for substring in substrings:
        n_chars = len(substring)

        for i in range(n_chars - 1):
            start_point = special_points[substring[i]]
            end_point = special_points[substring[i + 1]]
            start_point_str = format_coordinate(start_point)
            end_point_str = format_coordinate(end_point)
            line = f"    recipr    {start_point_str}    {end_point_str}    40\n"
            text += line

    text = f"$kpoints \n    kptlines {nlines}\n" + text
    return text

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

st.title('Band Structure Path')

structure = None

cif_file = st.file_uploader("Upload CIF file")
cif_contents = st.text_area("Or paste CIF contents here", height=300)

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
        visualize_structure(primitive_structure, "viz1.html")
        st.success("Converted to Primitive Structure! Using primitive structure from now on.")
        display_structure_info(primitive_structure)
        atoms = AseAtomsAdaptor.get_atoms(primitive_structure)
    else:
        st.warning("Using Conventional Structure. May result in Band Folding")
        visualize_structure(structure, "viz1.html")
        display_structure_info(structure)
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
    bandstructure_input = generate_turbomole_text(nlines, bandpath_str, special_points)
    # lattice_info_input = generate_lattice_text(structure)
    turbomole_text = st.text_area("`control` file text", value=bandstructure_input, height=200)



