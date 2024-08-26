import streamlit as st
import numpy as np
from io import StringIO
from ase import Atoms
from ase.build import make_supercell
from ase.io import read, write
from abtem.structures import orthogonalize_cell, is_cell_orthogonal
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
import py3Dmol
import pandas as pd
import streamlit.components.v1 as components

st.set_page_config(page_title='CIF Orthogonalization Tool', layout='wide', page_icon="⚛️")

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

def visualize_structure(structure, html_file_name='viz.html'):
    spin = st.checkbox('Spin', value=False, key='key' + html_file_name)
    view = py3Dmol.view(width=500, height=400)
    cif_for_visualization = structure.to(fmt="cif")
    view.addModel(cif_for_visualization, 'cif')
    view.setStyle({'sphere': {'colorscheme': 'Jmol', 'scale': 0.3},
                   'stick': {'colorscheme': 'Jmol', 'radius': 0.2}})
    view.addUnitCell()
    view.zoomTo()
    view.spin(spin)
    view.setClickable({'clickable': 'true'})
    view.enableContextMenu({'contextMenuEnabled': 'true'})
    view.show()
    view.render()
    t = view.js()
    f = open(html_file_name, 'w')
    f.write(t.startjs)
    f.write(t.endjs)
    f.close()

    HtmlFile = open(html_file_name, 'r', encoding='utf-8')
    source_code = HtmlFile.read()
    components.html(source_code, height=400, width=500)
    HtmlFile.close()

def display_structure_info(structure, atoms):
    st.subheader("Structure Information")
    st.write("Formula: ", structure.composition.reduced_formula)

    a, b, c = atoms.cell.lengths()
    alpha, beta, gamma = atoms.cell.angles()

    data = {
        "Lattice Parameters": [a, b, c, alpha, beta, gamma]
    }
    df_latt_params = pd.DataFrame(data, index=["a", "b", "c", "alpha", "beta", "gamma"])
    st.write("Lattice Parameters (Angstroms):")
    st.table(df_latt_params.style.format('{:.8f}'))

    lattice_vectors = atoms.cell[:]
    df_vectors = pd.DataFrame(lattice_vectors, columns=["X", "Y", "Z"], index=["a", "b", "c"])
    st.write("Lattice Vectors (Angstroms):")
    st.table(df_vectors.style.format('{:.8f}'))

    st.write("Atomic Coordinates (Fractional):")
    atomic_coords = []
    for site in structure.sites:
        atomic_coords.append([site.species_string] + list(site.frac_coords))
    df_coords = pd.DataFrame(atomic_coords, columns=["Element", "X", "Y", "Z"])
    st.table(df_coords.style.format('{:.8f}'))

def orthogonalize_cif(cif_content):
    # Read the structure from CIF content
    atoms = read(StringIO(cif_content), format='cif')
    
    # Check if the cell is already orthogonal
    if is_cell_orthogonal(atoms):
        st.warning("The cell is already orthogonal.")
        return atoms, None
    
    # Perform the conversion to an orthogonal cell
    orthogonal_atoms, strain = orthogonalize_cell(atoms, max_repetitions=5, return_transform=True)
    
    return orthogonal_atoms, strain

st.title('CIF Orthogonalization Tool')

cif_input = st.radio("Choose input method:", ('Paste CIF contents', 'Upload CIF file'))

if cif_input == 'Paste CIF contents':
    cif_contents = st.text_area("Paste CIF contents here", height=300)
else:
    cif_file = st.file_uploader("Upload CIF file", type=['cif'])
    if cif_file is not None:
        cif_contents = cif_file.getvalue().decode("utf-8")
    else:
        cif_contents = None

if cif_contents:
    st.subheader("Original Structure")
    original_atoms = read(StringIO(cif_contents), format='cif')
    original_structure = AseAtomsAdaptor.get_structure(original_atoms)
    
    col1, col2 = st.columns(2)
    with col1:
        visualize_structure(original_structure, "viz_original.html")
    with col2:
        display_structure_info(original_structure, original_atoms)
    
    if st.button("Orthogonalize Cell"):
        orthogonal_atoms, strain = orthogonalize_cif(cif_contents)
        
        if strain is not None:
            st.subheader("Orthogonalized Structure")
            orthogonal_structure = AseAtomsAdaptor.get_structure(orthogonal_atoms)
            
            col1, col2 = st.columns(2)
            with col1:
                visualize_structure(orthogonal_structure, "viz_orthogonal.html")
            with col2:
                display_structure_info(orthogonal_structure, orthogonal_atoms)
            
            st.write("Strain:")
            st.write(strain)
            
            # Save orthogonalized CIF to a StringIO object
            output = StringIO()
            write(output, orthogonal_atoms, format='cif')
            orthogonal_cif = output.getvalue()
            
            st.download_button(
                label="Download Orthogonalized CIF",
                data=orthogonal_cif,
                file_name="orthogonalized_structure.cif",
                mime="chemical/x-cif"
            )
        else:
            st.info("The cell is already orthogonal. No changes were made.")
