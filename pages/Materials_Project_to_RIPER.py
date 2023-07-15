import streamlit as st
from mp_api.client import MPRester
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import py3Dmol
import pandas as pd
import streamlit.components.v1 as components

# Set page config
st.set_page_config(page_title='Materials Project ➡️ RIPER', layout='wide', page_icon="⚛️",
menu_items={
         'About': "A web app to help you with DFT related calculations using the RIPER module of [TURBOMOLE](https://www.turbomole.org/)"
     })

# Sidebar stuff
st.sidebar.write('# About')
st.sidebar.write('### Made By [Manas Sharma](https://manas.bragitoff.com)')
st.sidebar.write('### In the group of [Prof. Dr. Marek Sierka](https://cmsg.uni-jena.de)')
st.sidebar.write('### *Powered by*')
st.sidebar.write('* [Py3Dmol](https://3dmol.csb.pitt.edu/) for Chemical System Visualizations')
st.sidebar.write('* [Streamlit](https://streamlit.io/) for making of the Web App')
st.sidebar.write('* [PyMatgen](https://pymatgen.org/) for Periodic Structure Representations')
st.sidebar.write('* [PubChempy](https://pypi.org/project/PubChemPy/1.0/) for Accessing the PubChem Database')
st.sidebar.write('* [MP-API](https://pypi.org/project/mp-api/) for Accessing the Materials Project Database')

# Set your Materials Project API key
api_key = st.secrets["MP_API"]

# Function to convert atomic coordinates to Bohr units
def convert_to_bohr(structure):
    coords = [(site.coords[0], site.coords[1], site.coords[2], site.species_string) for site in structure.sites]
    return [(x * 1.8897259886, y * 1.8897259886, z * 1.8897259886, element.lower()) for x, y, z, element in coords]

# Function to generate coordinate text
def generate_coord_text(coords_bohr):
    coord_text = "$coord\n"
    for coord in coords_bohr:
        coord_text += f"    {coord[0]:.8f}   {coord[1]:.8f}   {coord[2]:.8f}    {coord[3]}\n"
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
    spin = st.checkbox('Spin', value = False, key='key'+html_file_name)
    view = py3Dmol.view(width=500, height=400)
    cif_for_visualization = structure.to(fmt="cif")
    view.addModel(cif_for_visualization, 'cif')
    # view.setStyle({'stick': {'radius': 0.2}})
    view.setStyle({'sphere':{'colorscheme':'Jmol','scale':0.3},
                    'stick':{'colorscheme':'Jmol', 'radius':0.2}})
    view.addUnitCell()
    view.zoomTo()
    view.spin(spin)
    view.setClickable({'clickable':'true'});
    view.enableContextMenu({'contextMenuEnabled':'true'})
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
    components.html(source_code, height = 300, width=900)
    HtmlFile.close()
    

# Function to display structure information
def display_structure_info(structure):
    st.subheader("Structure Information")
    st.write("Formula: ", structure.composition.reduced_formula)

    #Display lattice parameters
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
        if coord_type=='Cartesian':
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

# Create an instance of MPRester
mpr = MPRester(api_key)

# Streamlit app
docs = None
st.write('# Materials Project ➡️ RIPER')
st.write("#### Get atomic coordinates and cell parameters for RIPER (TURBOMOLE) from Materials Project Database")

# Search for materials
formula = st.text_input("Enter formula:", placeholder='NaCl')
if not formula=="":
    with st.spinner("Searching..."):
        docs = mpr.summary.search(formula=[formula], fields=["structure", "band_gap", "material_id", "is_stable", "is_metal", "symmetry", "formula_pretty"])

        if len(docs) > 0:
            st.success(f"Matching materials found: {len(docs)}")

            # Display materials as a table
            table_data = []
            for doc in docs:
                table_data.append({
                    "Material ID": doc.material_id,
                    "Band Gap": doc.band_gap,
                    "Crystal System": str(doc.symmetry.crystal_system),
                    "Symbol": str(doc.symmetry.symbol),
                    "Group": str(doc.symmetry.point_group),
                    "Formula": doc.formula_pretty,
                    "Is Stable": doc.is_stable,
                    "Is Metal": doc.is_metal
                })
            st.table(table_data)
        else:
            st.warning("No materials found for the given formula.")

if docs is not None:

    # Select a material
    selected_material = st.selectbox("Select a material:", [doc.material_id for doc in docs])
    selected_doc = next((doc for doc in docs if doc.material_id == selected_material), None)
    structure = selected_doc.structure
    # Get conventional structure
    analyzer = SpacegroupAnalyzer(structure)
    conventional_structure = analyzer.get_conventional_standard_structure()
    # Get primitive structure
    primitive_structure = analyzer.get_primitive_standard_structure()

    # Choose between primitive or conventional
    selected_structure_type = st.selectbox("Unit cell type:", ['Primitve Cell','Conventional Unit Cell'])

    # Display structure information
    if selected_structure_type=='Primitive Cell':
        display_structure_info(primitive_structure)
    else:
        display_structure_info(conventional_structure)

    # Visualize the structure
    if selected_structure_type=='Primitive Cell':
        visualize_structure(primitive_structure, "viz1.html")
    else:
        visualize_structure(conventional_structure, "viz1.html")

    # Download CIF files
    st.subheader("Download CIF Files")
    col1, col2 = st.columns(2)
    convert_to_cif(primitive_structure, "primitive_unit_cell.cif")
    col1.download_button('Download Primitive Unit Cell CIF', "primitive_unit_cell.cif", file_name='primitive_unit_cell.cif', key='primitive_cif_button')
    st.success("Primitive Unit Cell CIF downloaded!")
    convert_to_cif(conventional_structure, "conventional_unit_cell.cif")
    col2.download_button('Download Conventional Unit Cell CIF', "conventional_unit_cell.cif", file_name='conventional_unit_cell.cif', key='conventional_cif_button')
    st.success("Conventional Unit Cell CIF downloaded!")

    # Get TURBOMOLE (RIPER) Coord file and Control file contents
    st.subheader("RIPER Files")
    # Convert the atomic coordinates to Bohr units
    if selected_structure_type=='Primitive Cell':
        coords_bohr = convert_to_bohr(primitive_structure)
    else:
        coords_bohr = convert_to_bohr(conventional_structure)

    # Generate the coordinate text
    coords_text = generate_coord_text(coords_bohr)
    # Generate the lattice parameter text
    if selected_structure_type=='Primitive Cell':
        lattice_text = generate_lattice_text(primitive_structure)
    else:
        lattice_text = generate_lattice_text(conventional_structure)

    # Create two columns for text boxes
    col1, col2 = st.columns(2)

    # Display the coordinate text in the first column
    with col1:
        st.text_area("Coord file contents (Cartesian coordinates in Bohr)", value=coords_text, height=300, key='coords_text')
        st.download_button('Download coord file', coords_text, file_name='coord', key='control_text')

    # Display the lattice parameters text in the second column
    with col2:
        st.text_area("Add the following to your control file", value=lattice_text, height=300)


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
    if selected_structure_type=='Primitive Cell':
        supercell_structure = primitive_structure.copy()
        supercell_structure.make_supercell([int(nx), int(ny), int(nz)])
    else:
        supercell_structure = conventional_structure.copy()
        supercell_structure.make_supercell([int(nx), int(ny), int(nz)])
    visualize_structure(supercell_structure, 'viz2.html')

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
        st.text_area("Coord file contents (Cartesian coordinates in Bohr)", value=coords_text_super, height=300, key='supercell_text_coord')
        st.download_button('Download coord file', coords_text_super, file_name='coord')

    # Display the lattice parameters text in the second column
    with col2:
        st.text_area("Add the following to your control file", value=lattice_text_super, height=300, key='supercell_text_control')

