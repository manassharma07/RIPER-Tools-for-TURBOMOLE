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
from rdkit import Chem
from rdkit.Chem import AllChem
from pymatgen.core.structure import Molecule
from mace.calculators import mace_mp

# Sidebar stuff
st.sidebar.write('# About')
st.sidebar.write(' Originally Made By [Manas Sharma](https://manas.bragitoff.com)')
st.sidebar.write(' In the group of [Prof. Dr. Marek Sierka](https://cmsg.uni-jena.de)')
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



def optimize_geometry_rdkit(smiles: str) -> Molecule:
    """
    Optimize geometry using RDKit's UFF implementation for a molecule given as SMILES.
    :param smiles: The SMILES string of the molecule.
    :return: Optimized geometry as a Pymatgen Molecule.
    """
    try:
        # Generate RDKit molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        # Add hydrogens
        mol = Chem.AddHs(mol)

        # Embed 3D coordinates
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        if result != 0:
            raise RuntimeError("3D embedding failed. Check the SMILES or molecular structure.")

        # Optimize geometry using UFF
        result = AllChem.UFFOptimizeMolecule(mol)
        if result != 0:
            # raise RuntimeError("UFF optimization failed.")
            st.warning('Optimization not converged after maxIterations!')

        # Extract optimized geometry
        conf = mol.GetConformer()
        optimized_coords = [
            [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
            for i in range(mol.GetNumAtoms())
        ]
        optimized_species = [atom.GetSymbol() for atom in mol.GetAtoms()]

        # Convert to Pymatgen Molecule
        optimized_molecule = Molecule(optimized_species, optimized_coords)
        return optimized_molecule

    except Exception as e:
        raise RuntimeError(f"Optimization failed: {e}")


# Function to format floating-point numbers with alignment
def format_number(num, width=10, precision=5):
    # Handles positive/negative numbers while maintaining alignment
    return f"{num: {width}.{precision}f}"

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
        xyz_text += f"{atom_symbol} {format_number(x, precision=8)} {format_number(y, precision=8)} {format_number(z, precision=8)}\n"

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
    components.html(source_code, height = 300, width=500)
    HtmlFile.close()
    

def search_pubchem(query):
    compounds = pcp.get_compounds(query, namespace='name', as_dataframe=False)
    return compounds


def get_molecule(cid):
    xyz_str = generate_xyz_coordinates(cid)
    return Molecule.from_str(xyz_str, fmt='xyz'), xyz_str


def format_xyz(molecule):
    xyz = molecule.to(fmt="xyz")
    return xyz


def format_coord(molecule):
    coord = ['$coord']
    for site in molecule:
        coord.append(
            f"   {format_number(site.x* 1.88972612456506, precision=8)}   {format_number(site.y* 1.88972612456506, precision=8)}   {format_number(site.z* 1.88972612456506, precision=8)}       {site.specie.symbol.lower()}"
        )
    coord.append('$end')
    return "\n".join(coord)

@st.cache_resource
def get_mace_mp():
    return mace_mp(model="https://github.com/ACEsuit/mace-mp/releases/download/mace_omat_0/mace-omat-0-medium.model", device="cpu", default_dtype="float32")
# device = "cpu" #"cuda"
# Initialize the MACE-MP calculator
# mace_mp_calc = mace_mp(model="small", device=device, default_dtype="float32")

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
    selected_molecule, xyz_str = get_molecule(selected_cid)
    selected_smiles = molecule_df.loc[molecule_df["CID"] == selected_cid, "Isomeric SMILES"].values[0]

    st.subheader("3D Atomic Coordinates")
    # Create a dataframe with atomic symbols and atomic coordinates
    st.dataframe(pd.DataFrame({"Atomic Symbol": selected_molecule.species, "X": selected_molecule.cart_coords[:, 0], 
                               "Y": selected_molecule.cart_coords[:, 1], "Z": selected_molecule.cart_coords[:, 2]}))


    opt_geom = st.checkbox(label= 'Optimize Geometry via ML (MACE MP Foundation) Model (Beta - Does not work well yet)', value=False)
    if opt_geom:
        ase_atoms = AseAtomsAdaptor().get_atoms(selected_molecule)
        # Set up the calculator for energy and forces using ML (or other forcefield)
        ase_atoms.calc = get_mace_mp()
        # Set up the optimizer (BFGS in this example)
        optimizer = BFGS(ase_atoms)

        # Run the optimization
        optimizer.run(fmax=0.05, steps=20)  # Adjust fmax value as needed

        # Display the energies and forces at each step
        for i, atoms in enumerate(optimizer.traj):
            energy = atoms.get_potential_energy()
            forces = atoms.get_forces()
            st.write(f"Step {i + 1}: Energy = {energy:.4f} eV, Forces = {forces}")

        # Visualization of the optimized structure
        st.subheader("Optimized Geometry")
        visualize_structure(AseAtomsAdaptor().get_molecule(ase_atoms), html_file_name="optimized_viz.html")

        st.write('Total Energy (eV) from [MACE ML model trained on OMAT24 from Meta](https://github.com/ACEsuit/mace-mp/releases/tag/mace_omat_0)')
        st.write(optimizer.atoms.get_potential_energy())
        st.write('Forces (eV/Angs) from [MACE ML model trained on OMAT24 from Meta](https://github.com/ACEsuit/mace-mp/releases/tag/mace_omat_0)')
        st.write(optimizer.atoms.get_forces())

        # Display and download optimized coordinates
        optimized_xyz = format_xyz(AseAtomsAdaptor().get_molecule(ase_atoms))
        optimized_coord = format_coord(AseAtomsAdaptor().get_molecule(ase_atoms))

        col1, col2 = st.columns(2)
        col1.subheader("Optimized XYZ Format")
        col1.code(optimized_xyz)
        col2.subheader("Optimized Turbomole Coord Format")
        col2.code(optimized_coord)

        col1.download_button(
            "Download Optimized XYZ",
            data=optimized_xyz,
            file_name="optimized_molecule_rdkit.xyz",
            mime="chemical/x-xyz",
        )

        col2.download_button(
            "Download Optimized Coord",
            data=optimized_coord,
            file_name="optimized_coord_rdkit",
            mime="text/plain",
        )


        # Get the optimized structure as Pymatgen structure
        selected_molecule = AseAtomsAdaptor().get_molecule(ase_atoms)
    # Visualization
    visualize_structure(selected_molecule)


    # XYZ and Coord files
    col1, col2 = st.columns(2)
    col1.subheader("XYZ Format")
    # col1.code(format_xyz(selected_molecule))
    col1.code(xyz_str)

    col2.subheader("Turbomole Coord Format")
    col2.code(format_coord(selected_molecule))

    col1.download_button(
        "Download XYZ",
        data=xyz_str,
        file_name="molecule.xyz",
        mime="chemical/x-xyz"
    )

    col2.download_button(
        "Download Coord",
        data=format_coord(selected_molecule),
        file_name="coord",
        mime="text/plain"
    )

    if st.button("Optimize Geometry with RDKit UFF"):
        with st.spinner("Optimizing geometry using RDKit UFF..."):
            try:
                optimized_molecule = optimize_geometry_rdkit(selected_smiles)
                st.success("Geometry optimization with RDKit UFF completed!")

                # Visualization of the optimized structure
                st.subheader("Optimized Geometry")
                visualize_structure(optimized_molecule, html_file_name="optimized_viz.html")

                # Display and download optimized coordinates
                optimized_xyz = format_xyz(optimized_molecule)
                optimized_coord = format_coord(optimized_molecule)

                col1, col2 = st.columns(2)
                col1.subheader("Optimized XYZ Format")
                col1.code(optimized_xyz)
                col2.subheader("Optimized Turbomole Coord Format")
                col2.code(optimized_coord)

                col1.download_button(
                    "Download Optimized XYZ",
                    data=optimized_xyz,
                    file_name="optimized_molecule_rdkit.xyz",
                    mime="chemical/x-xyz",
                )

                col2.download_button(
                    "Download Optimized Coord",
                    data=optimized_coord,
                    file_name="optimized_coord_rdkit",
                    mime="text/plain",
                )
            except Exception as e:
                st.error(f"Optimization failed: {str(e)}")
