import streamlit as st
import subprocess
import os
import tempfile
from pymatgen.core import Structure, Lattice
from pymatgen.io.cif import CifWriter
import py3Dmol
import numpy as np

# ---- Helper functions for TURBOMOLE (RIPER) output ----

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

# ---- End helper functions ----

st.title("Molecular Packing with Packmol")

st.markdown("""
### Instructions
1. **Cell Definition**: Choose how to define the unit cell:
   - **Upload a CIF file**: Use this if you want to add molecules on a substrate/monolayer or add solvent molecules to an existing structure.
   - **Specify cell parameters manually**: Use this if you want to pack molecules into an empty cell with known dimensions.
   - **Determine cell from density**: Use this if you know the desired density and number of molecules. The cell parameters will be computed automatically (cubic cell).
2. **Packing Range (Fractional)**: Optionally restrict where molecules are packed along each axis using fractional coordinate ranges. For example, if your CIF has a monolayer at 0.5 along the c-axis and you only want molecules between 0.5 and 0.7 in fractional coordinates along c, set the c-range to (0.5, 0.7).
3. **Upload one or more XYZ files** of the molecules you want to pack.
4. **Set the number of molecules** for each uploaded molecule type.
5. **Set the tolerance** (minimum distance between molecules in Å).
6. Click **Run Packmol** to generate the packed structure.
7. The resulting structure is shown as a 3D visualization and available for download as a CIF file and TURBOMOLE RIPER files.
""")

# ---- Feature 1: Cell definition method ----
cell_method = st.radio(
    "How would you like to define the unit cell?",
    ["Upload a CIF file", "Specify cell parameters manually", "Determine cell from density"],
    index=0
)

cif_structure = None
a, b, c = 10.0, 10.0, 10.0
alpha, beta, gamma = 90.0, 90.0, 90.0

if cell_method == "Upload a CIF file":
    cif_file = st.file_uploader("Upload CIF file", type=["cif"])
    if cif_file:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".cif") as tmp:
            tmp.write(cif_file.read())
            tmp_cif_path = tmp.name
        cif_structure = Structure.from_file(tmp_cif_path)
        a, b, c = cif_structure.lattice.abc
        alpha, beta, gamma = cif_structure.lattice.angles
        st.write(f"Loaded CIF with lattice parameters: a={a:.4f}, b={b:.4f}, c={c:.4f}, α={alpha:.2f}, β={beta:.2f}, γ={gamma:.2f}")

elif cell_method == "Specify cell parameters manually":
    col_a, col_b, col_c = st.columns(3)
    with col_a:
        a = st.number_input("a (Å)", value=10.0, min_value=0.1, format="%.4f")
    with col_b:
        b = st.number_input("b (Å)", value=10.0, min_value=0.1, format="%.4f")
    with col_c:
        c = st.number_input("c (Å)", value=10.0, min_value=0.1, format="%.4f")
    col_al, col_be, col_ga = st.columns(3)
    with col_al:
        alpha = st.number_input("α (°)", value=90.0, min_value=1.0, max_value=179.0, format="%.2f")
    with col_be:
        beta = st.number_input("β (°)", value=90.0, min_value=1.0, max_value=179.0, format="%.2f")
    with col_ga:
        gamma = st.number_input("γ (°)", value=90.0, min_value=1.0, max_value=179.0, format="%.2f")
    cif_structure = None

elif cell_method == "Determine cell from density":
    st.markdown("Provide the desired density, molecular weight, and number of molecules. A **cubic cell** will be generated.")
    density = st.number_input("Density (g/cm³)", value=1.0, min_value=0.01, format="%.4f")
    mol_weight = st.number_input("Molecular weight of one molecule (g/mol)", value=18.015, min_value=0.1, format="%.4f")
    num_mol_density = st.number_input("Number of molecules", value=10, min_value=1, step=1)

    # Calculate volume: V = (n * M) / (rho * N_A)
    # n = number of molecules, M = molecular weight (g/mol), rho = density (g/cm^3), N_A = Avogadro's number
    N_A = 6.02214076e23
    volume_cm3 = (num_mol_density * mol_weight) / (density * N_A)  # in cm^3
    volume_A3 = volume_cm3 * 1e24  # convert cm^3 to Å^3
    a = volume_A3 ** (1.0 / 3.0)
    b, c = a, a
    alpha, beta, gamma = 90.0, 90.0, 90.0
    st.write(f"Computed cubic cell: a = b = c = {a:.4f} Å, Volume = {volume_A3:.2f} ų")
    cif_structure = None

# Build the lattice from whatever method was chosen
lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

# ---- Feature 2: Packing range in fractional coordinates ----
st.subheader("Packing Range (Fractional Coordinates)")
st.markdown("Restrict where molecules are placed along each axis. Default is the full cell (0.0 to 1.0).")

col_fx, col_fy, col_fz = st.columns(3)
with col_fx:
    frac_x_min = st.number_input("a min (frac)", value=0.0, min_value=0.0, max_value=1.0, format="%.4f", key="fxmin")
    frac_x_max = st.number_input("a max (frac)", value=1.0, min_value=0.0, max_value=1.0, format="%.4f", key="fxmax")
with col_fy:
    frac_y_min = st.number_input("b min (frac)", value=0.0, min_value=0.0, max_value=1.0, format="%.4f", key="fymin")
    frac_y_max = st.number_input("b max (frac)", value=1.0, min_value=0.0, max_value=1.0, format="%.4f", key="fymax")
with col_fz:
    frac_z_min = st.number_input("c min (frac)", value=0.0, min_value=0.0, max_value=1.0, format="%.4f", key="fzmin")
    frac_z_max = st.number_input("c max (frac)", value=1.0, min_value=0.0, max_value=1.0, format="%.4f", key="fzmax")

# Convert fractional ranges to Cartesian ranges for Packmol
# For a general cell, we convert the 8 corners of the fractional box to Cartesian
# and use the axis-aligned bounding box
frac_corners = [
    [fx, fy, fz]
    for fx in [frac_x_min, frac_x_max]
    for fy in [frac_y_min, frac_y_max]
    for fz in [frac_z_min, frac_z_max]
]
cart_corners = np.array([lattice.get_cartesian_coords(fc) for fc in frac_corners])
pack_min = cart_corners.min(axis=0)  # [x_min, y_min, z_min]
pack_max = cart_corners.max(axis=0)  # [x_max, y_max, z_max]

st.write(f"Packing region (Cartesian, Å): x=[{pack_min[0]:.4f}, {pack_max[0]:.4f}], y=[{pack_min[1]:.4f}, {pack_max[1]:.4f}], z=[{pack_min[2]:.4f}, {pack_max[2]:.4f}]")

# ---- Molecule uploads ----
xyz_files = st.file_uploader("Upload XYZ file(s) of molecules", type=["xyz"], accept_multiple_files=True)

molecule_data = []
if xyz_files:
    for i, xyz_file in enumerate(xyz_files):
        st.write(f"**Molecule {i+1}: {xyz_file.name}**")
        num_molecules = st.number_input(f"Number of '{xyz_file.name}' molecules", value=10, min_value=1, step=1, key=f"num_mol_{i}")
        molecule_data.append((xyz_file, num_molecules))

tolerance = st.number_input("Tolerance (minimum distance between molecules, Å)", value=2.0, min_value=0.5, format="%.2f")

if st.button("Run Packmol"):
    if not molecule_data:
        st.error("Please upload at least one XYZ file.")
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Write molecule XYZ files to temp directory
            mol_paths = []
            for i, (xyz_file, num_mol) in enumerate(molecule_data):
                mol_path = os.path.join(tmpdir, f"molecule_{i}.xyz")
                with open(mol_path, "w") as f:
                    f.write(xyz_file.getvalue().decode("utf-8"))
                mol_paths.append(mol_path)

            output_path = os.path.join(tmpdir, "packed_output.xyz")

            # Build Packmol input
            packmol_input = f"tolerance {tolerance}\nfiletype xyz\noutput {output_path}\n\n"

            # If we have a CIF structure, fix those atoms first
            if cif_structure is not None:
                fixed_xyz_path = os.path.join(tmpdir, "fixed_structure.xyz")
                # Write the CIF structure atoms as a fixed XYZ
                with open(fixed_xyz_path, "w") as f:
                    f.write(f"{len(cif_structure)}\n")
                    f.write("Fixed atoms from CIF\n")
                    for site in cif_structure:
                        f.write(f"{site.species_string} {site.coords[0]:.6f} {site.coords[1]:.6f} {site.coords[2]:.6f}\n")
                packmol_input += f"structure {fixed_xyz_path}\n"
                packmol_input += f"  number 1\n"
                packmol_input += f"  fixed 0. 0. 0. 0. 0. 0.\n"
                packmol_input += f"end structure\n\n"

            # Add each molecule type with the packing region
            for i, (mol_path, (_, num_mol)) in enumerate(zip(mol_paths, molecule_data)):
                packmol_input += f"structure {mol_path}\n"
                packmol_input += f"  number {num_mol}\n"
                packmol_input += f"  inside box {pack_min[0]:.6f} {pack_min[1]:.6f} {pack_min[2]:.6f} {pack_max[0]:.6f} {pack_max[1]:.6f} {pack_max[2]:.6f}\n"
                packmol_input += f"end structure\n\n"

            # Write Packmol input file
            input_path = os.path.join(tmpdir, "packmol.inp")
            with open(input_path, "w") as f:
                f.write(packmol_input)

            st.text("Packmol Input:")
            st.code(packmol_input)

            # Run Packmol
            try:
                result = subprocess.run(
                    ["packmol"],
                    input=packmol_input,
                    capture_output=True,
                    text=True,
                    cwd=tmpdir,
                    timeout=120
                )
                st.text("Packmol Output:")
                st.code(result.stdout)
                if result.returncode != 0:
                    st.error("Packmol failed. See output above.")
                    st.code(result.stderr)
                else:
                    # Read output XYZ
                    if os.path.exists(output_path):
                        with open(output_path, "r") as f:
                            xyz_content = f.read()

                        # Parse XYZ to pymatgen Structure with the lattice
                        lines = xyz_content.strip().split("\n")
                        n_atoms = int(lines[0])
                        species = []
                        coords = []
                        for line in lines[2:2+n_atoms]:
                            parts = line.split()
                            species.append(parts[0])
                            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

                        packed_structure = Structure(
                            lattice,
                            species,
                            coords,
                            coords_are_cartesian=True
                        )

                        # Show 3D visualization
                        st.subheader("Packed Structure Visualization")
                        cif_writer = CifWriter(packed_structure)
                        cif_string = str(cif_writer)

                        view = py3Dmol.view(width=600, height=400)
                        view.addModel(cif_string, "cif")
                        view.setStyle({"stick": {}, "sphere": {"radius": 0.3}})
                        view.addUnitCell()
                        view.zoomTo()
                        st.components.v1.html(view._make_html(), height=420, width=620)

                        # Download CIF
                        st.download_button("Download packed structure as CIF", cif_string, file_name="packed_structure.cif")

                        # Download XYZ
                        st.download_button("Download packed structure as XYZ", xyz_content, file_name="packed_structure.xyz")

                        # ---- TURBOMOLE (RIPER) output ----
                        st.subheader("RIPER Files")
                        coords_bohr = convert_to_bohr(packed_structure)
                        coords_text = generate_coord_text(coords_bohr)
                        lattice_text = generate_lattice_text(packed_structure)

                        col1, col2 = st.columns(2)
                        with col1:
                            st.text_area("Coord file contents (Cartesian coordinates in Bohr)", value=coords_text, height=300)
                            st.download_button('Download coord file', coords_text, file_name='coord')
                        with col2:
                            st.text_area("Add the following to your control file", value=lattice_text, height=300)

                    else:
                        st.error("Output file not found.")
            except FileNotFoundError:
                st.error("Packmol is not installed or not found in PATH. Please install Packmol.")
            except subprocess.TimeoutExpired:
                st.error("Packmol timed out. Try reducing the number of molecules or increasing the cell size.")
