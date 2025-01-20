import streamlit as st
import io
from ase import Atoms
from ase.io import read, write
from ase.build import add_adsorbate
from pymatgen.io.ase import AseAtomsAdaptor
import py3Dmol
import pandas as pd
import streamlit.components.v1 as components
from pymatgen.core import Structure, Element, Molecule, Lattice

# Set page config
st.set_page_config(page_title='Pack Molecules in Cell', layout='wide', page_icon="⚛️",
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

# Function to format floating-point numbers with alignment
def format_number(num, width=15, precision=8):
    # Handles positive/negative numbers while maintaining alignment
    # Adjusting the width based on the sign to ensure alignment
    return f"{num:>{width}.{precision}f}"

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

@st.fragment
def download_packed_struture(packed_structure):
    cif_output = "total_structure.cif"
    write(cif_output, packed_structure, format='cif')
    with open(cif_output, "rb") as f:
        st.download_button(
            label="Download Final Structure (CIF)",
            data=f,
            file_name="total_structure.cif",
            mime="chemical/x-cif"
        )

# Function to visualize the structure using py3Dmol
# @st.fragment
def visualize_structure(structure, column, html_file_name='viz.html'):
    with column:
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
@st.fragment
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
        coord_type = st.selectbox('Coordinate type', ['Cartesian', 'Fractional/Crystal'], key='selectbox'+structure.composition.reduced_formula)
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

st.title("Place Adsorbate (molecule) on a Surface (periodic)")
with st.expander('How to Use?', expanded=False):
    st.write("1. The tool works as follows.")
    st.write("   - Provide a `CIF` file of the cell (surface) you want to place the molecule on.")
    st.write("   - Provide the `XYZ` file of the molecule to place on the surface.")
    st.write("2. For both files, you can either:")
    st.write("   - Paste the file contents directly in the text area.")
    st.write("   - Upload the file using the file uploader.")
    st.write("3. Adjust the x, y, z position of the adsorbate based on your preference.")

st.info('The tool assumes that the surface CIF provided has the vacuum along the z direction, that is the surface is in the xy plane.')
col1, col2 = st.columns(2)
# CIF input section
col1.header("Base Structure (CIF)")
cif_input_method = col1.radio("Choose input method for CIF", ["Upload File", "Paste Content"], key="cif_method", horizontal=True)

base_structure = None
if cif_input_method == "Paste Content":
    cif_content = col1.text_area(
        label="Enter the contents of the CIF file",
        value="",
        placeholder="Paste your CIF file contents here...",
        height=200,
        key="cif_text_area"
    )
    if cif_content.strip():
        try:
            base_structure = read(io.StringIO(cif_content), format="cif")
            col1.success("Base structure loaded from pasted content.")
        except Exception as e:
            col1.error(f"Error reading CIF content: {str(e)}")
else:
    cif_file = col1.file_uploader("Upload CIF file", type=["cif"], key="cif_uploader")
    if cif_file:
        try:
            base_structure = read(cif_file, format="cif")
            col1.success("Base structure loaded from file.")
        except Exception as e:
            col1.error(f"Error reading CIF file: {str(e)}")

# XYZ input section
col2.header("Molecule Structure (XYZ)")
xyz_input_method = col2.radio("Choose input method for XYZ", ["Upload File", "Paste Content"], key="xyz_method", horizontal=True)

molecule = None
if xyz_input_method == "Paste Content":
    xyz_content = col2.text_area(
        label="Enter the contents of the XYZ file",
        value="",
        placeholder="Paste your XYZ file contents here...",
        height=200,
        key="xyz_text_area"
    )
    if xyz_content.strip():
        try:
            molecule = read(io.StringIO(xyz_content), format="xyz")
            col2.success("Molecule structure loaded from pasted content.")
        except Exception as e:
            col2.error(f"Error reading XYZ content: {str(e)}")
else:
    xyz_file = col2.file_uploader("Upload XYZ file", type=["xyz"], key="xyz_uploader")
    if xyz_file:
        try:
            molecule = read(io.StringIO(xyz_file.read().decode("utf-8")), format="xyz")
            col2.success("Molecule structure loaded from file.")
        except Exception as e:
            col2.error(f"Error reading XYZ file: {str(e)}")

# Parameters and adsorbate placement section
if base_structure is not None and molecule is not None:
    col1.write(f"Number of atoms in base structure: {len(base_structure)}")
    col2.write(f"Number of atoms in molecule: {len(molecule)}")
    # Get pymatgen structure for further processing if needed
    # base_structure_pymatgen = AseAtomsAdaptor().get_structure(base_structure)
    # molecule_pymatgen = AseAtomsAdaptor().get_molecule(molecule)
    # Translate molecule so its center of mass (COM) is at the origin
    com_mol = molecule.get_center_of_mass()
    # molecule.translate(-com_mol)
    st.write(com_mol)
    
    
    # Set up adsorbate parameters
    col2.subheader("Adjust Adsorbate Position")
    translate_x = col2.slider("Translate adsorbate along x (fractional coordinates)", min_value=0.0, max_value=1.0, step=0.01, value=0.5)
    translate_y = col2.slider("Translate adsorbate along y (fractional coordinates)", min_value=0.0, max_value=1.0, step=0.01, value=0.5)
    # translate_z = st.slider("Translate adsorbate along z (fractional coordinates)", min_value=0.0, max_value=1.0, step=0.01, value=0.5)

    # Rotation sliders
    rotate_x = col2.slider("Rotate molecule around x-axis (degrees)", min_value=0, max_value=360, step=1)
    rotate_y = col2.slider("Rotate molecule around y-axis (degrees)", min_value=0, max_value=360, step=1)
    rotate_z = col2.slider("Rotate molecule around z-axis (degrees)", min_value=0, max_value=360, step=1)

    # Apply rotations to molecule
    # molecule.rotate(rotate_x, 'x', center=(0, 0, 0))
    # molecule.rotate(rotate_y, 'y', center=(0, 0, 0))
    # molecule.rotate(rotate_z, 'z', center=(0, 0, 0))
    molecule.rotate(rotate_x, 'x', center=com_mol)
    molecule.rotate(rotate_y, 'y', center=com_mol)
    molecule.rotate(rotate_z, 'z', center=com_mol)
    
    # Convert fractional coordinates to Cartesian and apply translation
    adsorbate_position = (translate_x * base_structure.cell[0] +
                          translate_y * base_structure.cell[1])

    # Add adsorbate onto the surface at a specified height
    adsorbate_height = col2.slider("Adsorbate Height (Å)", min_value=-10.0, max_value=15.0, value=2.0, step=0.1)
    add_adsorbate(base_structure, molecule, adsorbate_height, position=adsorbate_position[:2]+com_mol[:2])
    packed_structure_pymatgen = AseAtomsAdaptor().get_structure(base_structure)
    
    col1.subheader("Structure Preview and Download")
    visualize_structure(packed_structure_pymatgen, col1, "viz1.html")
    # Display packed structure info
    col1.success("The molecule has been positioned on the surface. You can download the final structure as a CIF file.")
    # Download option
    download_packed_struture(base_structure)

    display_structure_info(packed_structure_pymatgen)

    # Get TURBOMOLE (RIPER) Coord file and Control file contents
    st.subheader("RIPER Files")
    # Convert the atomic coordinates to Bohr units
    coords_bohr = convert_to_bohr(packed_structure_pymatgen)

    # Generate the coordinate text
    coords_text = generate_coord_text(coords_bohr)

    # Generate the lattice parameter text
    if isinstance(packed_structure_pymatgen, Structure):
        lattice_text = generate_lattice_text(packed_structure_pymatgen)

    col1_, col2_ = st.columns(2)
    # Display the coordinate text in the first column
    with col1_:
        st.text_area("Coord file contents (Cartesian coordinates in Bohr)", value=coords_text, height=300,
                     key='coords_text')
        st.download_button('Download coord file', coords_text, file_name='coord', key='control_text')

    if isinstance(packed_structure_pymatgen, Structure):
        # Display the lattice parameters text in the second column
        with col2_:
            st.text_area("Add the following to your control file", value=lattice_text, height=300)
