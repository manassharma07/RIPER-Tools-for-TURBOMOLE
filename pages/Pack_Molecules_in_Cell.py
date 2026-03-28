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
    import random
    import string
    
    # Generate a random 6-character string
    random_suffix = ''.join(random.choices(string.ascii_letters + string.digits, k=6))

    with st.expander("Atomic Coordinates", expanded=False):
        coord_type = st.selectbox('Coordinate type', ['Cartesian', 'Fractional/Crystal'], key='selectbox'+structure.composition.reduced_formula+ '_' + random_suffix)
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


# Function to convert atomic coordinates to Bohr units
def convert_to_bohr(structure):
    coords = [(site.coords[0], site.coords[1], site.coords[2], site.species_string) for site in structure.sites]
    return [(x * 1.88972612456506, y * 1.88972612456506, z * 1.88972612456506, element.lower()) for x, y, z, element in coords]

# Function to generate coordinate text
def generate_coord_text(coords_bohr):
    coord_text = "$coord\n"
    for coord in coords_bohr:
        coord_text += f"    {coord[0]:.10f}   {coord[1]:.10f}   {coord[2]:.10f}    {coord[3]}\n"
    coord_text += "$end"
    return coord_text

# Function to generate lattice parameter text
def generate_lattice_text(structure):
    lattice_params = structure.lattice.abc
    angles = structure.lattice.angles
    lattice_text = "$cell angs\n"
    lattice_text += f"  {lattice_params[0]:.10f}   {lattice_params[1]:.10f}   {lattice_params[2]:.10f}   {angles[0]}   {angles[1]}   {angles[2]}\n"
    lattice_text += "$periodic 3\n"
    lattice_text += "$kpoints\n"
    lattice_text += "    nkpoints 1 1 1 # Gamma point calculation"
    return lattice_text


def random_rotation_matrix():
    """Generate a random rotation matrix using scipy Rotation."""
    return R.random().as_matrix()


def _get_all_image_positions(positions, cell_matrix):
    """Replicate positions into 27 periodic images (3x3x3) for minimum image overlap check."""
    shifts = np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1])).T.reshape(-1, 3)
    # shifts is (27, 3) in fractional space
    cart_shifts = shifts @ cell_matrix  # (27, 3) in Cartesian
    # positions is (N, 3); result is (27*N, 3)
    all_positions = (positions[np.newaxis, :, :] + cart_shifts[:, np.newaxis, :]).reshape(-1, 3)
    return all_positions


def has_overlap(base_positions, mol_positions, tolerance, cell_matrix):
    """Check if there's an overlap between base structure atoms and molecule atoms
    using KDTree with 27-image replication for correct periodic boundary handling."""
    if len(base_positions) == 0:
        return False

    # Replicate base positions into 27 images so that minimum-image distances are captured
    base_images = _get_all_image_positions(base_positions, cell_matrix)

    # Build tree on the (larger) replicated base; query with the intact molecule
    base_tree = cKDTree(base_images)
    distances, _ = base_tree.query(mol_positions, distance_upper_bound=tolerance)
    return np.any(distances < tolerance)


def pack_structure(base_structure, molecule, num_molecules, tolerance, frac_range=None):
    """Pack molecule into the base structure without overlaps.
    
    frac_range: optional dict with keys 'a', 'b', 'c', each being a tuple (min, max) 
                in fractional coordinates to constrain molecule placement.
    """
    # start_time = time.time()  # Start timing
    
    # Create a mutable copy of the base structure
    packed_structure = base_structure.copy()
    original_positions = molecule.get_positions()
    
    # Get molecule center for proper rotation
    mol_center = original_positions.mean(axis=0)
    cell_matrix = np.array(base_structure.get_cell())  # 3x3 matrix
    inv_cell = np.linalg.inv(cell_matrix)

    # Keep a running array of all base atom Cartesian positions (avoids repeated get_positions)
    all_base_positions = packed_structure.get_positions().copy()

    max_attempts = 500  # Limit to avoid infinite loops
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
                
                # 2. Apply random rotation (correct column-vector convention)
                rotation_matrix = random_rotation_matrix()
                rotated_positions = (rotation_matrix @ centered_positions.T).T
                
                # 3. Generate random displacement uniformly in fractional space, convert to Cartesian
                # If frac_range is specified, constrain the fractional coordinates
                if frac_range is not None:
                    frac_displacement = np.array([
                        np.random.uniform(frac_range['a'][0], frac_range['a'][1]),
                        np.random.uniform(frac_range['b'][0], frac_range['b'][1]),
                        np.random.uniform(frac_range['c'][0], frac_range['c'][1]),
                    ])
                else:
                    frac_displacement = np.random.rand(3)
                displacement = frac_displacement @ cell_matrix
                
                # 4. Apply displacement and set new positions
                final_positions = rotated_positions + displacement
                current_molecule.set_positions(final_positions)
                
                # Check for overlaps using the intact molecule positions (no per-atom PBC wrap)
                if not has_overlap(all_base_positions, final_positions, tolerance, cell_matrix):
                    packed_structure += current_molecule.copy()  # Add molecule to the packed structure
                    # Update the running base positions array
                    all_base_positions = np.vstack([all_base_positions, final_positions])
                    added = True
                    st.write(f"Molecule copy #{i+1} added successfully at position {displacement} after {attempt} attempts.")
            if added:
                print(f"Added the {i+1} th molecule at {attempt}th attempt")
            else:
                st.write(f"Failed to add the {i+1}th copy of the molecule after {max_attempts} attempts.")
                # st.write(f"Failed to add the {i+1}th copy of the molecule after {max_attempts} attempts. Stopping.")
                # break  # No point continuing; the cell is likely full
    
    return packed_structure

@st.fragment
def download_packed_struture(packed_structure):
    import tempfile
    import os

    # CIF
    with tempfile.NamedTemporaryFile(suffix='.cif', delete=False) as tmp:
        tmp_cif = tmp.name
    write(tmp_cif, packed_structure, format='cif')
    with open(tmp_cif, 'rb') as f:
        st.download_button(
            label="Download Packed Structure (CIF)",
            data=f,
            file_name="packed_structure.cif",
            mime="chemical/x-cif"
        )
    os.unlink(tmp_cif)

    # Extended XYZ
    with tempfile.NamedTemporaryFile(suffix='.xyz', delete=False) as tmp:
        tmp_extxyz = tmp.name
    write(tmp_extxyz, packed_structure, format='extxyz')
    with open(tmp_extxyz, 'rb') as f:
        st.download_button(
            label="Download Packed Structure (extXYZ)",
            data=f,
            file_name="packed_structure.xyz",
            mime="chemical/x-xyz"
        )
    os.unlink(tmp_extxyz)

    # POSCAR
    with tempfile.NamedTemporaryFile(suffix='.vasp', delete=False) as tmp:
        tmp_poscar = tmp.name
    write(tmp_poscar, packed_structure, format='vasp')
    with open(tmp_poscar, 'rb') as f:
        st.download_button(
            label="Download Packed Structure (POSCAR)",
            data=f,
            file_name="POSCAR",
            mime="text/plain"
        )
    os.unlink(tmp_poscar)


def estimate_cell_from_density(molecule, num_molecules, density):
    """Estimate a cubic cell parameter from target density and number of molecules.
    
    density: target density in g/cm^3
    Returns cell edge length in Angstroms for a cubic cell.
    """
    from ase.data import atomic_masses, atomic_numbers
    # Calculate molecular mass in g/mol
    mol_mass = sum(atomic_masses[atomic_numbers[sym]] for sym in molecule.get_chemical_symbols())
    # Total mass in grams
    total_mass_g = (num_molecules * mol_mass) / 6.02214076e23
    # Volume in cm^3
    volume_cm3 = total_mass_g / density
    # Volume in Angstrom^3 (1 cm = 1e8 Angstrom)
    volume_A3 = volume_cm3 * 1e24
    # Cubic cell edge
    cell_edge = volume_A3 ** (1.0 / 3.0)
    return cell_edge


# # Streamlit interface

st.title("Pack Molecules in a Cell")
st.write('1. The tool works as follows.')
st.write('   - Provide a `CIF` file of the cell (e.g., a substrate or monolayer) in which you want to pack molecules, **OR** specify the cell parameters manually, **OR** specify the target density and number of molecules to auto-determine a cubic cell.')
st.write('   - Provide the `XYZ` file of the molecule to pack the cell with.')
st.write('2. For both the cell and the molecule, you can either:')
st.write('   - Paste the file contents directly in the text area')
st.write('   - Upload the file using the file uploader')
st.write('3. Specify the number of molecules to be packed and the tolerance distance.')
st.write('4. Optionally, restrict the packing region using fractional coordinate ranges along a, b, and c.')

# Base structure input section
st.header("Base Structure")
cell_input_method = st.radio(
    "Choose how to define the cell",
    ["CIF File", "Manual Cell Parameters", "Auto from Density"],
    key="cell_method"
)

base_structure = None

if cell_input_method == "CIF File":
    # CIF input section
    st.subheader("CIF Input")
    cif_input_method = st.radio("Choose input method for CIF", ["Paste Content", "Upload File"], key="cif_method")

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

elif cell_input_method == "Manual Cell Parameters":
    st.subheader("Cell Parameters")
    col_a, col_b, col_c = st.columns(3)
    with col_a:
        cell_a = st.number_input("a (Å)", min_value=0.1, value=10.0, step=0.1, key="cell_a")
    with col_b:
        cell_b = st.number_input("b (Å)", min_value=0.1, value=10.0, step=0.1, key="cell_b")
    with col_c:
        cell_c = st.number_input("c (Å)", min_value=0.1, value=10.0, step=0.1, key="cell_c")
    col_alpha, col_beta, col_gamma = st.columns(3)
    with col_alpha:
        cell_alpha = st.number_input("α (°)", min_value=1.0, max_value=179.0, value=90.0, step=1.0, key="cell_alpha")
    with col_beta:
        cell_beta = st.number_input("β (°)", min_value=1.0, max_value=179.0, value=90.0, step=1.0, key="cell_beta")
    with col_gamma:
        cell_gamma = st.number_input("γ (°)", min_value=1.0, max_value=179.0, value=90.0, step=1.0, key="cell_gamma")
    
    # Create an empty ASE Atoms object with the specified cell
    from ase.cell import Cell
    cell = Cell.fromcellpar([cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma])
    base_structure = Atoms(cell=cell, pbc=True)
    st.success(f"Empty cell created: a={cell_a}, b={cell_b}, c={cell_c}, α={cell_alpha}, β={cell_beta}, γ={cell_gamma}")

elif cell_input_method == "Auto from Density":
    st.subheader("Density-based Cell")
    st.write("Specify the target density and the number of molecules. A cubic cell will be determined automatically.")
    st.info("You need to provide the molecule (XYZ) below first, then come back and click 'Generate Cell'.")
    auto_density = st.number_input("Target density (g/cm³)", min_value=0.01, value=1.0, step=0.01, key="auto_density")
    auto_num_molecules = st.number_input("Number of molecules (for cell estimation)", min_value=1, value=10, step=1, key="auto_num_mol")
    # base_structure will be set after molecule is loaded (see below)

if base_structure is not None and cell_input_method != "Auto from Density":
    base_structure_pymatgen = AseAtomsAdaptor().get_structure(base_structure)
    if len(base_structure) > 0:
        display_structure_info(base_structure_pymatgen)
        visualize_structure(base_structure_pymatgen, "viz1.html")
    else:
        st.write("Empty cell (no atoms). Lattice parameters:")
        a, b, c = base_structure_pymatgen.lattice.abc
        alpha, beta, gamma = base_structure_pymatgen.lattice.angles
        st.write(f"a={a:.4f} Å, b={b:.4f} Å, c={c:.4f} Å, α={alpha:.1f}°, β={beta:.1f}°, γ={gamma:.1f}°")

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

# Handle "Auto from Density" cell generation after molecule is loaded
if cell_input_method == "Auto from Density" and molecule is not None:
    if st.button("Generate Cell from Density"):
        cell_edge = estimate_cell_from_density(molecule, auto_num_molecules, auto_density)
        base_structure = Atoms(cell=[cell_edge, cell_edge, cell_edge], pbc=True)
        st.success(f"Cubic cell generated: a = b = c = {cell_edge:.4f} Å (density ≈ {auto_density} g/cm³ for {auto_num_molecules} molecules)")
        st.session_state['auto_base_structure'] = base_structure

# Retrieve auto-generated base structure from session state if applicable
if cell_input_method == "Auto from Density" and base_structure is None:
    if 'auto_base_structure' in st.session_state:
        base_structure = st.session_state['auto_base_structure']

# Parameters and packing section
if base_structure is not None and molecule is not None:
    st.write(f"Number of atoms in base structure: {len(base_structure)}")
    st.write(f"Number of atoms in molecule: {len(molecule)}")
    if len(base_structure)<=800 and len(molecule)<=40:
        st.header("Packing Parameters")
        num_molecules = st.slider(
            "Number of molecules to add",
            min_value=1,
            max_value=200,
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

        # Fractional coordinate range for packing region
        st.subheader("Packing Region (Fractional Coordinates)")
        use_frac_range = st.checkbox("Restrict packing to a specific region?", value=False, key="use_frac_range")
        frac_range = None
        if use_frac_range:
            col_fa, col_fb, col_fc = st.columns(3)
            with col_fa:
                frac_a_min = st.number_input("a min", min_value=0.0, max_value=1.0, value=0.0, step=0.01, key="frac_a_min")
                frac_a_max = st.number_input("a max", min_value=0.0, max_value=1.0, value=1.0, step=0.01, key="frac_a_max")
            with col_fb:
                frac_b_min = st.number_input("b min", min_value=0.0, max_value=1.0, value=0.0, step=0.01, key="frac_b_min")
                frac_b_max = st.number_input("b max", min_value=0.0, max_value=1.0, value=1.0, step=0.01, key="frac_b_max")
            with col_fc:
                frac_c_min = st.number_input("c min", min_value=0.0, max_value=1.0, value=0.0, step=0.01, key="frac_c_min")
                frac_c_max = st.number_input("c max", min_value=0.0, max_value=1.0, value=1.0, step=0.01, key="frac_c_max")
            frac_range = {
                'a': (frac_a_min, frac_a_max),
                'b': (frac_b_min, frac_b_max),
                'c': (frac_c_min, frac_c_max),
            }

        isWrap = st.toggle("Wrap structure?", True)

        if st.button("Pack Structure"):
            overall_start_time = time.time()
            
            with st.spinner("Packing molecules..."):
                try:
                    packed_structure = pack_structure(base_structure, molecule, num_molecules, tolerance, frac_range=frac_range)
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

                    # TURBOMOLE RIPER output
                    st.subheader("RIPER Files")
                    coords_bohr = convert_to_bohr(packed_structure_pymatgen)
                    coords_text = generate_coord_text(coords_bohr)
                    lattice_text = generate_lattice_text(packed_structure_pymatgen)

                    col1, col2 = st.columns(2)
                    with col1:
                        st.text_area("Coord file contents (Cartesian coordinates in Bohr)", value=coords_text, height=300, key="riper_coord")
                        st.download_button('Download coord file', coords_text, file_name='coord', key="dl_coord")
                    with col2:
                        st.text_area("Add the following to your control file", value=lattice_text, height=300, key="riper_control")

                except Exception as e:
                    st.error(f"Error during packing: {str(e)}")
    else:
        st.info('The web app can only allow using base structure with less than 800 atoms and molecule with less than 40 atoms. This is because of limited computational capacity of the server. For no cocnstraints download the source code from GitHub and run the code locally.')
else:
    st.info("Please provide both the base structure (cell) and molecule structure (XYZ) to continue.")
