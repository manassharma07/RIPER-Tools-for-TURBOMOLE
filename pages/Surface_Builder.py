import streamlit as st
import io
from ase import Atoms
from ase.io import read, write
from pymatgen.io.ase import AseAtomsAdaptor
import py3Dmol
import pandas as pd
import streamlit.components.v1 as components
from pymatgen.core import Structure, Element, Molecule, Lattice
from ase.build import surface, make_supercell

# Set page config
st.set_page_config(page_title='Build surface/slab', layout='wide', page_icon="⚛️",
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
        coord_type = st.selectbox('Coordinate type', ['Cartesian', 'Fractional/Crystal'], key='selectbox'+structure.composition.formula)
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

# Function to check if the structure is a bulk or slab
def is_bulk(structure):
    lattice = structure.lattice.matrix
    lengths = structure.lattice.abc
    # If any lattice vector is significantly larger, it's likely a slab with vacuum
    return all(length > 30 for length in lengths)

st.title("Build Slab/Surface")
# with st.expander('How to Use?', expanded=False):
#     st.write("1. The tool works as follows.")
#     st.write("   - Provide a `CIF` file of the cell (surface) you want to place the molecule on.")
#     st.write("   - Provide the `XYZ` file of the molecule to place on the surface.")
#     st.write("2. For both files, you can either:")
#     st.write("   - Paste the file contents directly in the text area.")
#     st.write("   - Upload the file using the file uploader.")
#     st.write("3. Adjust the x, y, z position of the adsorbate based on your preference.")

st.info('The tool assumes that the CIF provided corresponds to the conventional unit cell of the bulk material.')
# CIF input section
st.header("Bulk Structure (CIF)")
cif_input_method = st.radio("Choose input method for CIF", ["Upload File", "Paste Content"], key="cif_method", horizontal=True)

bulk_structure = None
if cif_input_method == "Paste Content":
    cif_content = st.text_area(
        label="Enter the contents of the CIF file",
        value="",
        placeholder="Paste your CIF file contents here...",
        height=200,
        key="cif_text_area"
    )
    if cif_content.strip():
        try:
            bulk_structure = read(io.StringIO(cif_content), format="cif")
            st.success("Base structure loaded from pasted content.")
        except Exception as e:
            st.error(f"Error reading CIF content: {str(e)}")
else:
    cif_file = st.file_uploader("Upload CIF file", type=["cif"], key="cif_uploader")
    if cif_file:
        try:
            bulk_structure = read(cif_file, format="cif")
            st.success("Base structure loaded from file.")
        except Exception as e:
            st.error(f"Error reading CIF file: {str(e)}")

if bulk_structure is None:
    st.error("Please provide a CIF file for the bulk structure.")
    st.stop()
pymatgen_structure = AseAtomsAdaptor.get_structure(bulk_structure)

st.subheader("Structure Info")
st.write(f"**Formula:** {pymatgen_structure.composition.reduced_formula}")
display_structure_info(pymatgen_structure)

if is_bulk(pymatgen_structure):
        st.success("The uploaded structure is a bulk material.")
else:
    st.warning("The uploaded structure is a slab.")

# Miller indices input
miller_indices = st.text_input("Specify Miller Indices (e.g., 0 0 1):")

# Vacuum size input
vacuum_size = st.slider("Vacuum size (Å):", min_value=0.0, max_value=50.0, value=15.0)

# Number of layers input
layers = st.slider("Number of layers:", min_value=1, max_value=15, value=5)

if st.button("Generate Surface Slab"):
    if miller_indices:
        miller = tuple(map(int, miller_indices.split()))
        slab = surface(bulk_structure, miller, layers, vacuum=vacuum_size)
        st.success("Surface slab generated successfully.")
        slab_pymatgen = AseAtomsAdaptor.get_structure(slab)

        # Visualization
        with st.expander("Visualize Slab Structure"):
            view = py3Dmol.view(width=800, height=400)
            view.addModel(slab_pymatgen.to(fmt="cif"), "cif")
            view.setStyle({'stick': {}})
            view.zoomTo()
            view.show()
            components.html(view._make_html(), height=400)

        # Download option
        # cif_output = "surface_slab.cif"
        # write(cif_output, slab, format="cif")
        # with open(cif_output, "rb") as f:
        #     st.download_button("Download Slab CIF", f, file_name="surface_slab.cif")
        # Download option
        download_packed_struture(slab)
        display_structure_info(slab_pymatgen)

        # Get TURBOMOLE (RIPER) Coord file and Control file contents
        st.subheader("RIPER Files")
        # Convert the atomic coordinates to Bohr units
        coords_bohr = convert_to_bohr(slab_pymatgen)

        # Generate the coordinate text
        coords_text = generate_coord_text(coords_bohr)

        # Generate the lattice parameter text
        if isinstance(slab_pymatgen, Structure):
            lattice_text = generate_lattice_text(slab_pymatgen)

        col1_, col2_ = st.columns(2)
        # Display the coordinate text in the first column
        with col1_:
            st.text_area("Coord file contents (Cartesian coordinates in Bohr)", value=coords_text, height=300,
                        key='coords_text')
            st.download_button('Download coord file', coords_text, file_name='coord', key='control_text')

        if isinstance(slab_pymatgen, Structure):
            # Display the lattice parameters text in the second column
            with col2_:
                st.text_area("Add the following to your control file", value=lattice_text, height=300)

    else:
        st.error("Please specify Miller indices.")

    