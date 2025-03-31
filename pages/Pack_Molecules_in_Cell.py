import streamlit as st
from ase import Atoms
from ase.io import read, write
import numpy as np
import io
import time  # For timing debug statements
from pymatgen.io.ase import AseAtomsAdaptor
import py3Dmol
import pandas as pd
import streamlit.components.v1 as components

from scipy.spatial.transform import Rotation as R
from scipy.spatial import cKDTree

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

# def has_overlap(base_structure, molecule, tolerance):
#     """Check if there's an overlap by calculating distances between the base structure and molecule."""
#     combined = base_structure + molecule
#     distances = combined.get_all_distances(mic=True)
#     np.fill_diagonal(distances, np.inf)
#     return np.min(distances) < tolerance

# def random_rotation_matrix():
#     """Generate a random rotation matrix using quaternions."""
#     # Random rotation angles
#     angles = np.random.rand(3) * 2 * np.pi
    
#     # Create rotation matrix using quaternions for better numerical stability
#     qx = np.array([[1, 0, 0],
#                    [0, np.cos(angles[0]), -np.sin(angles[0])],
#                    [0, np.sin(angles[0]), np.cos(angles[0])]])
    
#     qy = np.array([[np.cos(angles[1]), 0, np.sin(angles[1])],
#                    [0, 1, 0],
#                    [-np.sin(angles[1]), 0, np.cos(angles[1])]])
    
#     qz = np.array([[np.cos(angles[2]), -np.sin(angles[2]), 0],
#                    [np.sin(angles[2]), np.cos(angles[2]), 0],
#                    [0, 0, 1]])
    
#     return qz @ qy @ qx

# def has_overlap(base_structure, molecule, tolerance):
#     """Check if there's an overlap between base structure atoms and molecule atoms."""
#     # Get number of atoms in each structure
#     n_base = len(base_structure)
#     n_mol = len(molecule)
    
#     # Get all distances between base structure and molecule atoms only
#     combined = base_structure + molecule
#     distances = combined.get_all_distances(mic=True)
    
#     # Extract only the distances between base structure and molecule
#     # This is the submatrix of size (n_base x n_mol) from the distances matrix
#     relevant_distances = distances[:n_base, n_base:]
    
#     # Check if any distance is less than tolerance
#     return np.min(relevant_distances) < tolerance

def random_rotation_matrix():
    """Generate a random rotation matrix using scipy Rotation."""
    return R.random().as_matrix()

def has_overlap(base_structure, molecule, tolerance, cell):
    """Check if there's an overlap between base structure atoms and molecule atoms using KDTree with periodic boundary conditions."""
    base_positions = base_structure.get_positions()
    mol_positions = molecule.get_positions()
    
    if len(base_positions) == 0:
        return False
    
    # Use minimum image convention (mic) for periodic systems
    base_tree = cKDTree(base_positions, boxsize=cell)
    distances, _ = base_tree.query(mol_positions, distance_upper_bound=tolerance)
    return np.any(distances < tolerance)

def pack_structure(base_structure, molecule, num_molecules, tolerance):
    """Pack molecule into the base structure without overlaps."""
    # start_time = time.time()  # Start timing
    
    # Create a mutable copy of the base structure
    packed_structure = base_structure.copy()
    original_positions = molecule.get_positions()
    
    # Get molecule center for proper rotation
    mol_center = original_positions.mean(axis=0)
    cell = base_structure.get_cell().lengths()

    max_attempts = 50  # Limit to avoid infinite loops
    with st.expander("Packing...", expanded=False):
        # Loop to add molecules
        for i in range(num_molecules):
            attempt = 0
            added = False
            while not added and attempt < max_attempts:
                print(i, attempt)
                attempt += 1
                # Create a new copy of the molecule for this attempt
                current_molecule = molecule.copy()
                
                # 1. Center the molecule at origin
                centered_positions = original_positions - mol_center
                
                # 2. Apply random rotation
                rotation_matrix = random_rotation_matrix()
                rotated_positions = centered_positions @ rotation_matrix
                
                # 3. Generate random displacement within cell
                displacement = np.random.rand(3) * packed_structure.cell.lengths()
                
                # 4. Apply displacement and set new positions
                final_positions = rotated_positions + displacement
                current_molecule.set_positions(final_positions)
                
                # Check for overlaps
                if not has_overlap(packed_structure, current_molecule, tolerance, cell):
                    packed_structure += current_molecule.copy()  # Add molecule to the packed structure
                    added = True
                    st.write(f"Molecule copy #{i+1} added successfully at position {displacement} after {attempt} attempts.")
                # else:
                    # st.write(f"Attempt {attempt} for molecule {i+1} resulted in overlap; retrying.")
                # if attempt == max_attempts - 1:
                #     i = i - 1
            if added:
                print(f"Added the {i+1} th molecule at {attempt}th attempt")
            else:
                st.write(f"Failed to add the {i+1}th copy of the molecule after {max_attempts} attempts.")
                # Debug: Time for each molecule addition
                # st.write(f"Time to add molecule {i+1}: {time.time() - start_time:.2f} seconds")
                # start_time = time.time()  # Reset start time for next molecule
    
    return packed_structure

@st.fragment
def download_packed_struture(packed_structure):
    cif_output = "packed_structure.cif"
    write(cif_output, packed_structure, format='cif')
    with open(cif_output, "rb") as f:
        st.download_button(
            label="Download Packed Structure (CIF)",
            data=f,
            file_name="packed_structure.cif",
            mime="chemical/x-cif"
        )

# # Streamlit interface

st.title("Pack Molecules in a Cell")
st.write('1. The tool works as follows.')
st.write('   - Provide a `CIF` file of the cell in which you want to pack a number of some specific molecule.')
st.write('   - Provide the `XYZ` file of the molecule to pack the cell with.')
st.write('2. For both the cell and the molecule, you can either:')
st.write('   - Paste the file contents directly in the text area')
st.write('   - Upload the file using the file uploader')
st.write('3. Specify the number of molecules to be packed and the tolerance distance.')

# CIF input section
st.header("Base Structure (CIF)")
cif_input_method = st.radio("Choose input method for CIF", ["Paste Content", "Upload File"], key="cif_method")

base_structure = None
if cif_input_method == "Paste Content":
    cif_content = st.text_area(
        label='Enter the contents of the CIF file',
        value='',
        placeholder='Paste your CIF file contents here...',
        height=200,
        key='cif_text_area'
    )
    if cif_content.strip():
        try:
            base_structure = read(io.StringIO(cif_content), format='cif')
            st.success("Base structure loaded from pasted content.")
        except Exception as e:
            st.error(f"Error reading CIF content: {str(e)}")
else:
    cif_file = st.file_uploader("Upload CIF file", type=["cif"], key="cif_uploader")
    if cif_file:
        try:
            base_structure = read(cif_file, format='cif')
            st.success("Base structure loaded from file.")
        except Exception as e:
            st.error(f"Error reading CIF file: {str(e)}")
if base_structure is not None:
    base_structure_pymatgen = AseAtomsAdaptor().get_structure(base_structure)
    display_structure_info(base_structure_pymatgen)
    visualize_structure(base_structure_pymatgen, "viz1.html")
# XYZ input section
st.header("Molecule Structure (XYZ)")
xyz_input_method = st.radio("Choose input method for XYZ", ["Paste Content", "Upload File"], key="xyz_method")

molecule = None
if xyz_input_method == "Paste Content":
    xyz_content = st.text_area(
        label='Enter the contents of the XYZ file',
        value='',
        placeholder='Paste your XYZ file contents here...',
        height=200,
        key='xyz_text_area'
    )
    if xyz_content.strip():
        try:
            molecule = read(io.StringIO(xyz_content), format='xyz')
            st.success("Molecule structure loaded from pasted content.")
        except Exception as e:
            st.error(f"Error reading XYZ content: {str(e)}")
else:
    xyz_file = st.file_uploader("Upload XYZ file", type=["xyz"], key="xyz_uploader")
    if xyz_file:
        try:
            xyz_content = xyz_file.getvalue().decode("utf-8")
            molecule = read(io.StringIO(xyz_content), format='xyz')
            st.success("Molecule structure loaded from file.")
        except Exception as e:
            st.error(f"Error reading XYZ file: {str(e)}")
if molecule is not None:
    molecule_pymatgen = AseAtomsAdaptor().get_molecule(molecule)
    # display_structure_info(molecule)
    visualize_molecule(molecule_pymatgen, "viz2.html")

# Parameters and packing section
if base_structure is not None and molecule is not None:
    st.write(f"Number of atoms in base structure: {len(base_structure)}")
    st.write(f"Number of atoms in molecule: {len(molecule)}")
    if len(base_structure)<=800 and len(molecule)<=20:
        st.header("Packing Parameters")
        num_molecules = st.slider(
            "Number of molecules to add",
            min_value=1,
            max_value=100,
            value=1,
            step=1
        )
        tolerance = st.slider(
            "Tolerance distance between atoms (Å)",
            min_value=1.0,
            max_value=5.0,
            value=2.0,
            step=0.1
        )

        isWrap = st.toggle("Wrap structure?", True)

        if st.button("Pack Structure"):
            overall_start_time = time.time()
            
            with st.spinner("Packing molecules..."):
                try:
                    packed_structure = pack_structure(base_structure, molecule, num_molecules, tolerance)
                    st.success("Packed structure generated successfully!")
                    st.write(f"Total packing time: {time.time() - overall_start_time:.2f} seconds")
                    if packed_structure is not None:
                        
                        if isWrap:
                            packed_structure.wrap()
                        packed_structure_pymatgen = AseAtomsAdaptor().get_structure(packed_structure)
                        # packed_structure_pymatgen.apply_strain(1.0)
                        display_structure_info(packed_structure_pymatgen)
                        visualize_structure(packed_structure_pymatgen, "viz3.html")

                    # Download option
                    download_packed_struture(packed_structure)
                except Exception as e:
                    st.error(f"Error during packing: {str(e)}")
    else:
        st.info('The web app can only allow using base structure with less than 800 atoms and molecule with less than 20 atoms. This is because of limited computational capacity of the server. For no cocnstraints download the source code from GitHub and run the code locally.')
else:
    st.info("Please provide both the base structure (CIF) and molecule structure (XYZ) to continue.")
