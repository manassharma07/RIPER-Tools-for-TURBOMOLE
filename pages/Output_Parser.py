import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO
from pymatgen.core import Structure, Lattice, Site, Element
import py3Dmol
# from pymatgen.core.periodic_table import Element

def find_line_with_text(lines, text):
    for line in lines:
        if text in line:
            return line
    return None


def parse_energies(text):
    kinetic_energy = []
    coulomb_energy = []
    exchange_corr_energy = []
    total_energy = []
    
    lines = text.split('\n')
    for line in lines:
        if "KINETIC ENERGY" in line:
            kinetic_energy.append(float(line.split()[4]))
        elif "COULOMB ENERGY" in line:
            coulomb_energy.append(float(line.split()[4]))
        elif "EXCH. & CORR. ENERGY" in line:
            exchange_corr_energy.append(float(line.split()[6]))
        elif "TOTAL ENERGY" in line:
            total_energy.append(float(line.split()[4]))
    
    return kinetic_energy, coulomb_energy, exchange_corr_energy, total_energy

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

st.title("`RIPER` Output Parser")

latt_param_a = None
latt_param_b = None
latt_param_c = None
latt_param_alpha = None
latt_param_beta = None
latt_param_gamma = None

st.write('You can either paste the output file contents below or upload the source file')
contents = st.text_area(label='Enter the contents of the output file here', value='', placeholder='Put your text here',
                        height=400, key='input_text_area')
# Create a file uploader widget
file = st.file_uploader("or Upload the file")

if file is not None:
    # If a file is uploaded, read its contents
    # contents = file.read()
    # To read file as bytes:
    bytes_data = file.getvalue()

    # To convert to a string based IO:
    stringio = StringIO(file.getvalue().decode("utf-8"))

    # To read file as string:
    contents = stringio.read()

if contents != '':
    file_contents = contents #upload_file.read()
    energies = parse_energies(file_contents)
    
    st.subheader("Parsed Energies")
    data = {
        "SCF Iteration": list(range(1, len(energies[0]) + 1)),
        "Kinetic Energy": energies[0],
        "Coulomb Energy": energies[1],
        "Exchange Corr. Energy": energies[2],
        "Total Energy": energies[3]
    }
    df = pd.DataFrame(data)
    st.dataframe(df)

    st.subheader("Convergence (Energy vs SCF Iteration)")
    plt.figure(figsize=(10, 6))
    plt.plot(data["SCF Iteration"], data["Total Energy"], marker='o', linestyle='-', color='b', label='Total Energy')
    plt.xlabel("SCF Iteration")
    plt.ylabel("Energy")
    plt.title("Energy vs SCF Iteration")
    plt.legend()
    st.pyplot(plt)


    # Find periodicity and lattice parameters
    lines = file_contents.split('\n')
    cell_params_line = find_line_with_text(lines, "Cell parameters (au,deg.)")
    
    if cell_params_line is not None:
        st.success("The output file indicates a periodic DFT calculation and the structure can be visualized!")
        fourth_line = lines[lines.index(cell_params_line) + 4]
        num_elements = len(fourth_line.split())
        
        if num_elements == 6:
            periodicity = 3
        elif num_elements == 3:
            periodicity = 2
        elif num_elements == 1:
            periodicity = 1

        st.write('#### Periodicity: '+str(periodicity))

        if periodicity==1:
            latt_param_a = fourth_line.split()[0]
        if periodicity==2:
            latt_param_a = fourth_line.split()[0]
            latt_param_b = fourth_line.split()[1]
            latt_param_gamma = fourth_line.split()[2]
        if periodicity==3:
            latt_param_a = fourth_line.split()[0]
            latt_param_b = fourth_line.split()[1]
            latt_param_c = fourth_line.split()[2]
            latt_param_alpha = fourth_line.split()[3]
            latt_param_beta = fourth_line.split()[4]
            latt_param_gamma = fourth_line.split()[5]
            
        lattice_lines = []
        direct_space_line = find_line_with_text(lines, "Direct space cell vectors (au):")
        if direct_space_line is not None:
            lattice_lines = lines[lines.index(direct_space_line) + 1 : lines.index(direct_space_line) + periodicity + 1]
            

        # Find atomic coordinates
        fractional_coords_line = find_line_with_text(lines, "fractional coordinates")
        if fractional_coords_line is not None:
            atomic_coords_lines = lines[lines.index(fractional_coords_line) + 1:]
            atomic_coords = []
            
            for line in atomic_coords_lines:
                if line.strip() == "":
                    break  # Stop when an empty line is encountered
                parts = line.split()
                element = parts[0].capitalize()  # Assuming lowercase element symbols
                coords = list(map(float, parts[1:4]))
                atomic_coords.append((element, coords))

            # Create the lattice using the lattice vectors
            lattice_vectors = []
            for line in lattice_lines:
                lattice_vectors.append(list(map(float, line.split()[1:4])))
            lattice = Lattice(lattice_vectors)

            # Create the sites using atomic coordinates
            sites = []
            for element, coords in atomic_coords:
                sites.append({"species": element, "xyz": coords})

            # Create the Structure object
            structure = Structure(
                lattice=lattice,
                species=[site["species"] for site in sites],
                coords=[site["xyz"] for site in sites]
            )
            # Display structure information
            if isinstance(structure, Structure):  # type = Structure
                display_structure_info(structure)
                visualize_structure(structure, "viz1.html")

    else:
        st.error("Only output files of periodic DFT calculations can be visualized for now!")

