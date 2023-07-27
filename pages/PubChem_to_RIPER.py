import streamlit as st
import pubchempy as pcp
from pymatgen.core.structure import Molecule
import py3Dmol
import pandas as pd
import streamlit.components.v1 as components
from pymatgen.io.ase import AseAtomsAdaptor
from ase.calculators.emt import EMT
from ase.calculators.singlepoint import SinglePointCalculator
from ase.optimize import BFGS

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
st.sidebar.write('* [ASE](https://wiki.fysik.dtu.dk/ase/) for File Format Conversions')

# CID to XYZ
def generate_xyz_coordinates(cid):
    compound = pcp.Compound.from_cid(cid, record_type='3d')
    atoms = compound.atoms
    coords = [(atom.x, atom.y, atom.z) for atom in atoms]

    num_atoms = len(atoms)
    xyz_text = f"{num_atoms}\n{compound.cid}\n"

    for atom, coord in zip(atoms, coords):
        atom_symbol = atom.element
        x, y, z = coord
        xyz_text += f"{atom_symbol} {x} {y} {z}\n"

    return xyz_text


# Function to visualize the structure using py3Dmol
def visualize_structure(structure, html_file_name='viz.html'):
    spin = st.checkbox('Spin', value = False, key='key'+html_file_name)
    view = py3Dmol.view(width=500, height=400)
    xyz_for_visualization = structure.to(fmt="xyz")
    view.addModel(xyz_for_visualization, 'xyz')
    # view.setStyle({'stick': {'radius': 0.2}})
    view.setStyle({'sphere':{'colorscheme':'Jmol','scale':0.3},
                    'stick':{'colorscheme':'Jmol', 'radius':0.2}})
    # view.addUnitCell()
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
    

def search_pubchem(query):
    compounds = pcp.get_compounds(query, namespace='name', as_dataframe=False)
    return compounds


def get_molecule(cid):
    xyz_str = generate_xyz_coordinates(cid)
    return Molecule.from_str(xyz_str, fmt='xyz')


def format_xyz(molecule):
    xyz = molecule.to(fmt="xyz")
    return xyz


def format_coord(molecule):
    coord = ['$coord']
    for site in molecule:
        coord.append(
            f"   {site.x* 1.8897259886:.6f}   {site.y* 1.8897259886:.6f}   {site.z* 1.8897259886:.6f}       {site.specie.symbol.lower()}"
        )
    coord.append('$end')
    return "\n".join(coord)



st.title("PubChem ➡️ RIPER")

search_query = st.text_input("Enter molecule name or formula to search in the PubChem database", placeholder='Water / H2O')
molecule_df = None
compounds = None

if not search_query=="":
    compounds = search_pubchem(search_query)
    # st.dataframe(compounds) # if using the pubchempy dataframe

if compounds:
    molecule_df = pd.DataFrame(
        [(compound.cid, compound.iupac_name, compound.molecular_formula, compound.molecular_weight, compound.isomeric_smiles)
            for compound in compounds],
        columns=["CID", "Name", "Formula", "Weight", "Isomeric SMILES"]
    )

    st.success(f"{len(molecule_df)} molecule(s) found!")
    st.dataframe(molecule_df)

if compounds is not None:
    # selected_cid = st.selectbox("Select a molecule", [cid for cid in compounds.index]) # if using the pubchempy dataframe
    selected_cid = st.selectbox("Select a molecule", molecule_df["CID"])
    # st.write(generate_xyz_coordinates(selected_cid))
    selected_molecule = get_molecule(selected_cid)

    st.subheader("3D Atomic Coordinates")
    # Create a dataframe with atomic symbols and atomic coordinates
    st.dataframe(pd.DataFrame({"Atomic Symbol": selected_molecule.species, "X": selected_molecule.cart_coords[:, 0], 
                               "Y": selected_molecule.cart_coords[:, 1], "Z": selected_molecule.cart_coords[:, 2]}))


    # opt_geom = st.checkbox(label= 'Optimize Geometry via ASE EMT calculator (Beta - Does not work well yet)', value=False)
    # if opt_geom:
    #     ase_atoms = AseAtomsAdaptor().get_atoms(selected_molecule)
    #     # calc = SinglePointCalculator(ase_atoms, EMT())
    #     # Set up the calculator for energy and forces using EMT (or other forcefield)
    #     calc = EMT()
    #     ase_atoms.set_calculator(calc)
    #     # Set up the optimizer (BFGS in this example)
    #     optimizer = BFGS(ase_atoms)

    #     # Run the optimization
    #     optimizer.run(fmax=0.05, steps=10)  # Adjust fmax value as needed

    #     # Get the optimized structure as Pymatgen structure
    #     selected_molecule = AseAtomsAdaptor().get_molecule(ase_atoms)
    # Visualization
    visualize_structure(selected_molecule)


    # XYZ and Coord files
    col1, col2 = st.columns(2)
    col1.subheader("XYZ Format")
    col1.code(format_xyz(selected_molecule))

    col2.subheader("Turbomole Coord Format")
    col2.code(format_coord(selected_molecule))

    col1.download_button(
        "Download XYZ",
        data=format_xyz(selected_molecule),
        file_name="molecule.xyz",
        mime="chemical/x-xyz"
    )

    col2.download_button(
        "Download Coord",
        data=format_coord(selected_molecule),
        file_name="coord",
        mime="text/plain"
    )
