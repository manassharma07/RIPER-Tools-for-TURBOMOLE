import streamlit as st
import io
from ase import Atoms
from ase.io import read, write
from ase.build import add_adsorbate
# from pymatgen.io.ase import AseAtomsAdaptor
import tempfile

st.title("Place Adsorbate (molecule) on a Surface (periodic)")
st.write("1. The tool works as follows.")
st.write("   - Provide a `CIF` file of the cell (surface) you want to place the molecule on.")
st.write("   - Provide the `XYZ` file of the molecule to place on the surface.")
st.write("2. For both files, you can either:")
st.write("   - Paste the file contents directly in the text area.")
st.write("   - Upload the file using the file uploader.")
st.write("3. Adjust the x, y, z position of the adsorbate based on your preference.")

st.info('The tool assumes that the surface CIF provided has the vacuum along the z direction, that is the surface is in the xy plane.')

# CIF input section
st.header("Base Structure (CIF)")
cif_input_method = st.radio("Choose input method for CIF", ["Paste Content", "Upload File"], key="cif_method", horizontal=True)

base_structure = None
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
            base_structure = read(io.StringIO(cif_content), format="cif")
            st.success("Base structure loaded from pasted content.")
        except Exception as e:
            st.error(f"Error reading CIF content: {str(e)}")
else:
    cif_file = st.file_uploader("Upload CIF file", type=["cif"], key="cif_uploader")
    if cif_file:
        try:
            base_structure = read(cif_file, format="cif")
            st.success("Base structure loaded from file.")
        except Exception as e:
            st.error(f"Error reading CIF file: {str(e)}")

# XYZ input section
st.header("Molecule Structure (XYZ)")
xyz_input_method = st.radio("Choose input method for XYZ", ["Paste Content", "Upload File"], key="xyz_method")

molecule = None
if xyz_input_method == "Paste Content":
    xyz_content = st.text_area(
        label="Enter the contents of the XYZ file",
        value="",
        placeholder="Paste your XYZ file contents here...",
        height=200,
        key="xyz_text_area"
    )
    if xyz_content.strip():
        try:
            molecule = read(io.StringIO(xyz_content), format="xyz")
            st.success("Molecule structure loaded from pasted content.")
        except Exception as e:
            st.error(f"Error reading XYZ content: {str(e)}")
else:
    xyz_file = st.file_uploader("Upload XYZ file", type=["xyz"], key="xyz_uploader")
    if xyz_file:
        try:
            molecule = read(io.StringIO(xyz_file.read().decode("utf-8")), format="xyz")
            st.success("Molecule structure loaded from file.")
        except Exception as e:
            st.error(f"Error reading XYZ file: {str(e)}")

# Parameters and adsorbate placement section
if base_structure is not None and molecule is not None:
    # Get pymatgen structure for further processing if needed
    # base_structure_pymatgen = AseAtomsAdaptor().get_structure(base_structure)
    # molecule_pymatgen = AseAtomsAdaptor().get_molecule(molecule)
    # Translate molecule so its center of mass (COM) is at the origin
    molecule.translate(-molecule.get_center_of_mass())
    
    # Set up adsorbate parameters
    st.subheader("Adjust Adsorbate Position")
    translate_x = st.slider("Translate adsorbate along x (fractional coordinates)", min_value=0.0, max_value=1.0, step=0.01, value=0.5)
    translate_y = st.slider("Translate adsorbate along y (fractional coordinates)", min_value=0.0, max_value=1.0, step=0.01, value=0.5)
    translate_z = st.slider("Translate adsorbate along z (fractional coordinates)", min_value=0.0, max_value=1.0, step=0.01, value=0.5)

    # Rotation sliders
    rotate_x = st.slider("Rotate molecule around x-axis (degrees)", min_value=0, max_value=360, step=1)
    rotate_y = st.slider("Rotate molecule around y-axis (degrees)", min_value=0, max_value=360, step=1)
    rotate_z = st.slider("Rotate molecule around z-axis (degrees)", min_value=0, max_value=360, step=1)

    # Apply rotations to molecule
    molecule.rotate(rotate_x, 'x', center=(0, 0, 0))
    molecule.rotate(rotate_y, 'y', center=(0, 0, 0))
    molecule.rotate(rotate_z, 'z', center=(0, 0, 0))
    
    # Convert fractional coordinates to Cartesian and apply translation
    adsorbate_position = (translate_x * base_structure.cell[0] +
                          translate_y * base_structure.cell[1] +
                          translate_z * base_structure.cell[2])

    # Add adsorbate onto the surface at a specified height
    adsorbate_height = st.slider("Adsorbate Height (Å)", min_value=1.0, max_value=10.0, value=2.0, step=0.1)
    # add_adsorbate(base_structure, molecule, adsorbate_height, position=(translate_x * base_structure.cell[0][0:1] + translate_y * base_structure.cell[1][0:1]))
    add_adsorbate(base_structure, molecule, adsorbate_height, position=adsorbate_position[:2])
    
    # Display packed structure info
    st.subheader("Structure Preview and Download")
    st.write("The molecule has been positioned on the surface. You can download the final structure as a CIF file.")
    
    # Download button
    with tempfile.NamedTemporaryFile(delete=False, suffix=".cif") as tmp_file:
        write(tmp_file.name, base_structure, format="cif")
        tmp_file_path = tmp_file.name
    with open(tmp_file_path, "rb") as file:
        btn = st.download_button(
            label="Download Packed Structure as CIF",
            data=file,
            file_name="packed_structure.cif",
            mime="chemical/x-cif"
        )
