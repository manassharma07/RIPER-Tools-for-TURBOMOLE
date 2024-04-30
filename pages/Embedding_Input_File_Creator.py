import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
import subprocess
import sys
import time
from io import StringIO
import pandas as pd
import numpy as np
# from stmol import showmol

try:
    # from openbabel import OBMol, OBConversion
    import openbabel
except ModuleNotFoundError as e:
    subprocess.Popen([f'{sys.executable} -m pip install --global-option=build_ext --global-option="-I/usr/include/openbabel3" --global-option="-L/usr/lib/openbabel" openbabel'], shell=True)
    subprocess.Popen([f'{sys.executable} -m pip install --global-option=build_ext --global-option="-I/home/appuser/include/openbabel3" --global-option="-L/home/appuser/lib/openbabel" openbabel'], shell=True)
    subprocess.Popen([f'{sys.executable} -m pip install --global-option=build_ext --global-option="-I/home/appuser/usr/include/openbabel3" --global-option="-L/home/appuser/usr/lib/openbabel" openbabel'], shell=True)
    # wait for subprocess to install package before running your actual code below
    time.sleep(90)
    
import os
from openbabel import pybel

if os.path.exists('viz.html'):
    os.remove('viz.html')
if os.path.exists('viz1.html'):
    os.remove('viz1.html')

further_information = '''When you run the `riperembed.py` script, it reads in the input parameters from the file named 
'input' and the coordinates of the total system, from a file named `totalCoords`.\n
Next, it automatically creates the required subdirectories for the subsystems and runs `define` and `riper` automatically based on the given input file. \n
Let's have a look at the workflows of some specific casess.
##### Freeze-and-Thaw Workflow
If a Freeze-and-Thaw (FaT) calculation with `method=1,3,5 ` is requested then, it will create subdirectories `1`, `2`, `3`, ... corresponding 
to the `ith` FaT cycle.\n
1️⃣ For the first FaT cycle, the script will automatically run a regular KS-DFT calculation using `riper` on the isolated subsystem B (environment) 
within the directory `1` and store the output in a file called `output2`.\n
Note: The `control` and other TURBOMOLE related files are created by running `define` automatically. Then the script amends the `control` 
to add the embedding related keywords based on the input file.\n
Furthermore, auxiliary density coefficients of the subsystem B (`auxdcaoB`) and density matrix (`denssaoB`, `dmatspB`) are also saved on the disk. 
Some energies of the subsystem B are also saved in a file called `energyB`. These are in the following order XC energy, nuc-nuc energy and KEDF energies.\n
2️⃣ Then, the script will create a directory called `Cluster` (for subsystem A or active subsytem) and copy the necessary quantities, `denssaoB`, `dmatspB` 
and `auxdcaoB` to this directory and also `cd` to it. \n
Once again, `define` is run automatically and the necessary files are generated. Then the `control` file is modified with embedding related keywords for 
the active subsytem. \n
Subsequently, `riper` is launched and an embedding calculation is performed on the subsystem A in the presence of the isolated environment and the 
result is stored in `output1`.\n
3️⃣ Next, the script returns back to the `1` directory and creates another subdirectory called `Environment` for the subsystem B, and similar to before the 
required quantitities for embedding are copied from the `Cluster` directory and the embedding calculation is performed via `riper` on the subsystem B in 
the presence of subsystem B within the `Environment` directory and the result is stored in `output1`.\n
4️⃣ For the next FaT cycles, the steps 2️⃣-3️⃣ are repeated and the results are stored in directories `2`, `3`,... unless energies and densities are converged.
5️⃣ Finally, a regular KS-DFT calculation may be performed on the total system if requested in the input file, to estimate the embedding error.

'''

def COM_calculator(coords):
    return coords.mean(axis=0)

# Set page config
st.set_page_config(page_title="DFT based Embedding Input File Creator (for TURBOMOLE's `riper` module)", layout='wide', page_icon="🧊",
menu_items={
         'About': "### This online tool allows you to create an input file for DFT based embedding calculations using TURBOMOLE's riper module"
     })


# Sidebar stuff
st.sidebar.write('# About')
st.sidebar.write('### Made By [Manas Sharma](https://manas.bragitoff.com/)')
st.sidebar.write('### *Powered by*')
st.sidebar.write('* [Py3Dmol](https://pypi.org/project/py3Dmol/) for Visualization')
st.sidebar.write('* [Open Babel](http://openbabel.org/) for Format Conversion')
st.sidebar.write('## Brought to you by [TURBOMOLE](https://www.turbomole.org)')
st.sidebar.write('## Cite us:')
st.sidebar.write('[Manas Sharma and Marek Sierka, Journal of Chemical Theory and Computation **2022** 18 (11), 6892-6904](https://doi.org/10.1021/acs.jctc.2c00380)')
with st.sidebar.expander('Instructions'):
    st.write('1️⃣ Upload/Paste the contents of an XYZ file that contains the atomic coordinates of the total system.')
    st.write('2️⃣ Select the atoms belonging to subsystem A. These will automatically be placed in the beginning in the new coords file.')
    st.write('3️⃣ Download the totalCoords file.')
    st.write('4️⃣ Create and Download the input file using the GUI.')
    st.write('5️⃣ Put the two files (totalCoords adn input) in the same directory.')
    st.write('6️⃣ Run the riperembed.py script as: `nohup riperembed.py > output_embedding` &')

# Main app
st.write("## DFT based Embedding Input File Creator (for TURBOMOLE's riper module)")
st.write("This online tool allows you to create an input file for DFT based embedding calculations using TURBOMOLE's riper module")


# DATA for test systems
he_dimer_xyz = '''
2

He          1.00000        0.00000        0.00000
He          2.00000        0.00000        0.00000
'''
hf_dimer_xyz = '''
4

F          1.32374       -0.09023       -0.00001
H          1.74044        0.73339        0.00001
F         -1.45720        0.01926       -0.00001
H         -0.53931       -0.09466        0.00015
'''
h2o_dimer_xyz = '''
6

O          1.53175        0.00592       -0.12088
H          0.57597       -0.00525        0.02497
H          1.90625       -0.03756        0.76322
O         -1.39623       -0.00499        0.10677
H         -1.78937       -0.74228       -0.37101
H         -1.77704        0.77764       -0.30426
'''
nh3_dimer_xyz = '''
8

N          1.57523        0.00008       -0.04261
H          2.13111        0.81395       -0.28661
H          1.49645       -0.00294        0.97026
H          2.13172       -0.81189       -0.29145
N         -1.68825        0.00008        0.10485
H         -2.12640       -0.81268       -0.31731
H         -2.12744        0.81184       -0.31816
H         -0.71430        0.00054       -0.19241
'''
ch4_dimer_xyz = ''''''
benzene_fulvene_dimer_xyz = '''
24

C         -0.65914       -1.21034        3.98683
C          0.73798       -1.21034        4.02059
C         -1.35771       -0.00006        3.96990
C          1.43653       -0.00004        4.03741
C         -0.65915        1.21024        3.98685
C          0.73797        1.21024        4.02061
H         -1.20447       -2.15520        3.97369
H          1.28332       -2.15517        4.03382
H         -2.44839       -0.00006        3.94342
H          2.52722       -0.00004        4.06369
H         -1.20448        2.15509        3.97373
H          1.28330        2.15508        4.03386
C          0.64550       -0.00003        0.03038
C         -0.23458       -1.17916       -0.00274
C         -0.23446        1.17919       -0.00272
C         -1.51856       -0.73620       -0.05059
C         -1.51848        0.73637       -0.05058
C          1.99323       -0.00010        0.08182
H          0.11302       -2.20914        0.01010
H          0.11325        2.20913        0.01013
H         -2.41412       -1.35392       -0.08389
H         -2.41398        1.35418       -0.08387
H          2.56084        0.93137        0.10366
H          2.56074       -0.93163        0.10364
'''
ethane_xyz = '''
8

C          0.00000        0.00000        0.76510
H          0.00000       -1.02220        1.16660
H         -0.88530        0.51110        1.16660
H          0.88530        0.51110        1.16660
C          0.00000        0.00000       -0.76510
H          0.88530       -0.51110       -1.16660
H          0.00000        1.02220       -1.16660
H         -0.88530       -0.51110       -1.16660
'''

dict_name_to_xyz = {'HF dimer': hf_dimer_xyz,'He dimer': he_dimer_xyz, 'H2O dimer': h2o_dimer_xyz,'NH3 dimer': nh3_dimer_xyz,'Benzene-Fulvene': benzene_fulvene_dimer_xyz,'Ethane': ethane_xyz}

input_test_system = st.selectbox('Select a test system',
     ( 'HF dimer', 'He dimer', 'H2O dimer', 'NH3 dimer', 'Benzene-Fulvene', 'Ethane'))

selected_xyz_str = dict_name_to_xyz[input_test_system]

st.write('#### Alternatively you can provide the XYZ file contents of your own structure here')
st.write("If you don't have an XYZ file but some other format, then you can use this [online tool](https://crysx-file-converter.streamlitapp.com/) to convert it to XYZ")
input_text_area = st.empty()
input_geom_str = input_text_area.text_area(label='XYZ file of the given/selected system', value = selected_xyz_str, placeholder = 'Put your text here', height=250, key = 'input_text_area')
# Get rid of trailing empty lines
input_geom_str = input_geom_str.rstrip()
# Get rid of leading empty lines
input_geom_str = input_geom_str.lstrip()

uploaded_file = st.file_uploader("You can also choose a file on your system")
if uploaded_file is not None:
    # To read file as bytes:
    bytes_data = uploaded_file.getvalue()

    # To convert to a string based IO:
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))

    # To read file as string:
    string_data = stringio.read()
    selected_xyz_str = string_data
    input_geom_str = input_text_area.text_area(label='XYZ file of the given/selected system', value = selected_xyz_str, placeholder = 'Put your text here', height=250)
    
print(input_geom_str)
print(type(input_geom_str))
### Create a dataframe from the original XYZ file ###
# Separate into lines and remove the first two lines
inp_geom_str_splitlines = input_geom_str.splitlines()[2:]
INPUT_GEOM_DATA = StringIO("\n".join(inp_geom_str_splitlines))
df = pd.read_csv(INPUT_GEOM_DATA, delim_whitespace=True, names=['atom','x','y','z'])
# df.reindex(index=range(1, natoms_tot+1))
df.index += 1 
coords_Tot_np_arr = df[['x','y','z']].to_numpy()

st.write('#### Visualization')
### VISUALIZATION ####
style = st.selectbox('Visualization style',['ball-stick','line','cross','stick','sphere'])
col1, col2 = st.columns(2)
spin = col1.checkbox('Spin', value = False)
showLabels = col2.checkbox('Show Labels', value = True)
view = py3Dmol.view(width=500, height=300)
structure_for_visualization = ''
try:
    mol = pybel.readstring('xyz', input_geom_str)
    natoms_tot = len(mol.atoms)
    # mol.make3D()
    if style=='cartoon':
        structure_for_visualization = mol.write('pdb')
    else:
        structure_for_visualization = mol.write('xyz')
except Exception as e:
    print('There was a problem with the conversion', e)
if style=='cartoon':
    view.addModel(structure_for_visualization, 'pdb')
else:
    view.addModel(structure_for_visualization, 'xyz')
if style=='ball-stick': # my own custom style
    view.setStyle({'sphere':{'colorscheme':'Jmol','scale':0.3},
                       'stick':{'colorscheme':'Jmol', 'radius':0.2}})
else:
    view.setStyle({style:{'colorscheme':'Jmol'}})
# Label addition template
# view.addLabel('Aromatic', {'position': {'x':-6.89, 'y':0.75, 'z':0.35}, 
#             'backgroundColor': 'white', 'backgroundOpacity': 0.5,'fontSize':18,'fontColor':'black',
#                 'fontOpacity':1,'borderThickness':0.0,'inFront':'true','showBackground':'false'})
if showLabels:
    for atom in mol:
        view.addLabel(str(atom.idx), {'position': {'x':atom.coords[0], 'y':atom.coords[1], 'z':atom.coords[2]}, 
            'backgroundColor': 'white', 'backgroundOpacity': 0.5,'fontSize':18,'fontColor':'black',
                'fontOpacity':1,'borderThickness':0.0,'inFront':'true','showBackground':'false'})
# Draw Axis
originAxis_Offset = np.array([-2.0, -2.0, 1.0])
originAxis = originAxis_Offset + np.array([np.min(coords_Tot_np_arr[:,0]), np.min(coords_Tot_np_arr[:,1]), np.min(coords_Tot_np_arr[:,2])])
view.addArrow({"start": {"x":originAxis[0], "y":originAxis[1], "z":originAxis[2]}, "end": {"x":originAxis[0]+0.8, "y":originAxis[1], "z":originAxis[2]}, "radiusRadio": 0.2, "color":"red"})
view.addArrow({"start": {"x":originAxis[0], "y":originAxis[1], "z":originAxis[2]}, "end": {"x":originAxis[0], "y":originAxis[1]+0.8, "z":originAxis[2]}, "radiusRadio": 0.2, "color":"green"})
view.addArrow({"start": {"x":originAxis[0], "y":originAxis[1], "z":originAxis[2]}, "end": {"x":originAxis[0], "y":originAxis[1], "z":originAxis[2]+0.8}, "radiusRadio": 0.2, "color":"blue"})
# view.addLabel('x', {'position': {'x':originAxis[0]+0.8, 'y':originAxis[1] + 0.1, 'z':originAxis[2]}, 
#             'backgroundColor': 'white', 'backgroundOpacity': 0.5,'fontSize':15,'fontColor':'black',
#                 'fontOpacity':1,'borderThickness':0.0,'inFront':'true','showBackground':'false'})
# view.addLabel('y', {'position': {'x':originAxis[0], 'y':originAxis[1]+0.8 + 0.1, 'z':originAxis[2]}, 
#             'backgroundColor': 'white', 'backgroundOpacity': 0.5,'fontSize':15,'fontColor':'black',
#                 'fontOpacity':1,'borderThickness':0.0,'inFront':'true','showBackground':'false'})
# view.addLabel('z', {'position': {'x':originAxis[0], 'y':originAxis[1] + 0.1, 'z':originAxis[2]+0.8}, 
#             'backgroundColor': 'white', 'backgroundOpacity': 0.5,'fontSize':15,'fontColor':'black',
#                 'fontOpacity':1,'borderThickness':0.0,'inFront':'true','showBackground':'false'})
# view.addCylinder({"start": {"x":0.0, "y":2.0, "z":0.0}, "end": {"x":0.5, "y":2.0, "z":0.0}, "radius": 0.5, "color":"red"})
view.zoomTo()
view.spin(spin)
view.setClickable({'clickable':'true'})
view.enableContextMenu({'contextMenuEnabled':'true'})
view.show()
view.render()
view.resize()
# view.png()
t = view.js()
f = open('viz.html', 'w')
f.write(t.startjs)
f.write(t.endjs)
f.close()

col1, col2 = st.columns(2)
with col1:
    st.write('#### Visualization [Original]')
    HtmlFile = open("viz.html", 'r', encoding='utf-8')
    source_code = HtmlFile.read() 
    components.html(source_code, height = 300, width=900)
    HtmlFile.close()
    st.write('###### Axis labels')
    st.write('*x* : red, \n *y* :green, \n *z* :blue')

with col2:
    st.write('#### Atomic Positions ')
    st.dataframe(df, height=350)
    st.write('Center of Mass of the total system', COM_calculator(coords_Tot_np_arr))

# Select some rows using st.multiselect. This will break down when you have >1000 rows.
st.write('#### Choose the atom labels/indices that should assigend to subsystem A')
selected_indices = st.multiselect('Select rows:', df.index)
selected_rows_A = df.loc[selected_indices]
natoms_A = selected_rows_A.shape[0]
selected_rows_B = df.loc[~df.index.isin(selected_indices)]
# selected_rows_B.index -= natoms_A
natoms_B = selected_rows_B.shape[0]

coords_A_np_arr = selected_rows_A[['x','y','z']].to_numpy()
coords_B_np_arr = selected_rows_B[['x','y','z']].to_numpy()

if not natoms_A==0:
    with st.expander('More customizations', expanded=False):
        tab1, tab2 = st.tabs(['Translate', 'Rotate'])

        with tab1:
            subsystem_to_translate = st.selectbox('Choose a subsystem to translate', ['A','B'])
            col_translate1, col_translate2, col_translate3 = st.columns(3)
            translate_x = col_translate1.number_input('Translate along x', value=0.0, min_value=-10.0, max_value=10.0, step=0.1)
            translate_y = col_translate2.number_input('Translate along y', value=0.0, min_value=-10.0, max_value=10.0, step=0.1)
            translate_z = col_translate3.number_input('Translate along z', value=0.0, min_value=-10.0, max_value=10.0, step=0.1)
            translate_com = st.number_input('Translate along the line joining the COMs of the subsystems', value=0.0, min_value=-10.0, max_value=10.0, step=0.1)
            if subsystem_to_translate=='A':
                coords_A_np_arr[:,0] += translate_x
                coords_A_np_arr[:,1] += translate_y
                coords_A_np_arr[:,2] += translate_z
            if subsystem_to_translate=='B':
                coords_B_np_arr[:,0] += translate_x
                coords_B_np_arr[:,1] += translate_y
                coords_B_np_arr[:,2] += translate_z

        with tab2:
            subsystem_to_rotate = st.selectbox('Choose a subsystem to rotate', ['A','B'])
            col_translate1, col_translate2, col_translate3 = st.columns(3)
            translate_x = col_translate1.number_input('Rotate about x', value=0.0, min_value=0.0, max_value=360.0, step=1.0)
            translate_y = col_translate2.number_input('Rotate about y', value=0.0, min_value=0.0, max_value=360.0, step=1.0)
            translate_z = col_translate3.number_input('Rotate about z', value=0.0, min_value=0.0, max_value=360.0, step=1.0)
            translate_com = st.number_input('Rotate about the line joining the COMs of the subsystems', value=0.0, min_value=0.0, max_value=360.0, step=1.0)

    #Create an XYZ file of the modified/customized structure
    selected_rows_A.iloc[:,1:4] = coords_A_np_arr
    selected_rows_B.iloc[:,1:4] = coords_B_np_arr
    coords_Tot_np_arr_new = np.concatenate([coords_A_np_arr, coords_B_np_arr])
    modified_xyz = str(natoms_tot)+'\n Created for RIPER Embedding\n'
    modified_xyz = modified_xyz + selected_rows_A.to_string(header=False, index=False, float_format='{:.6f}'.format)
    modified_xyz = modified_xyz + '\n' + selected_rows_B.to_string(header=False, index=False, float_format='{:.6f}'.format)
    xyz_view = py3Dmol.view(width=500, height=300)
    xyz_view.addModel(modified_xyz, 'xyz')
    xyz_view.setStyle({'sphere':{'colorscheme':'Jmol','scale':0.3},
                        'stick':{'colorscheme':'Jmol', 'radius':0.2}})

    if showLabels:
        iat = 0
        for atom in mol:
            xyz_view.addLabel(str(atom.idx), {'position': {'x':coords_Tot_np_arr_new[iat,0], 'y':coords_Tot_np_arr_new[iat,1], 'z':coords_Tot_np_arr_new[iat,2]}, 
                'backgroundColor': 'white', 'backgroundOpacity': 0.5,'fontSize':18,'fontColor':'black',
                    'fontOpacity':1,'borderThickness':0.0,'inFront':'true','showBackground':'false'})
            iat = iat + 1
    # Draw Axis
    originAxis_Offset = np.array([-2.0, -2.0, 1.0])
    originAxis = originAxis_Offset + np.array([np.min(coords_Tot_np_arr_new[:,0]), np.min(coords_Tot_np_arr_new[:,1]), np.min(coords_Tot_np_arr_new[:,2])])
    xyz_view.addArrow({"start": {"x":originAxis[0], "y":originAxis[1], "z":originAxis[2]}, "end": {"x":originAxis[0]+0.8, "y":originAxis[1], "z":originAxis[2]}, "radiusRadio": 0.2, "color":"red"})
    xyz_view.addArrow({"start": {"x":originAxis[0], "y":originAxis[1], "z":originAxis[2]}, "end": {"x":originAxis[0], "y":originAxis[1]+0.8, "z":originAxis[2]}, "radiusRadio": 0.2, "color":"green"})
    xyz_view.addArrow({"start": {"x":originAxis[0], "y":originAxis[1], "z":originAxis[2]}, "end": {"x":originAxis[0], "y":originAxis[1], "z":originAxis[2]+0.8}, "radiusRadio": 0.2, "color":"blue"})

    xyz_view.zoomTo()
    xyz_view.setClickable({'clickable':'true'})
    xyz_view.enableContextMenu({'contextMenuEnabled':'true'})
    xyz_view.show()
    xyz_view.center()
    xyz_view.render()
    xyz_view.resize()

    t1 = xyz_view.js()
    f1 = open('viz1.html', 'w')
    f1.write(t1.startjs)
    f1.write(t1.endjs)
    f1.close()

    st.write('##### The following shows the reordered atom indices, such that the subsystem A atoms are in the beginning of the coordinate file')
    st.write("You can also translate or move subsystems using the 'More customizations' tab above, and the changes will reflect in the following visualization")
    HtmlFile1 = open("viz1.html", 'r', encoding='utf-8')
    source_code1 = HtmlFile1.read() 
    components.html(source_code1, height = 300, width=900)
    HtmlFile1.close()
    st.write('###### Axis labels')
    st.write('*x* : red, \n *y* :green, \n *z* :blue')

    com_A = COM_calculator(coords_A_np_arr)
    com_B = COM_calculator(coords_B_np_arr)
    dist_bw_COM_subsystems = np.linalg.norm(com_A - com_B)
    # The following is a 2d array that contains the euclidean distances between the atoms of the two subsystems
    dist_bw_atoms_subsystems = np.sqrt((coords_A_np_arr[:, 0, np.newaxis] - coords_B_np_arr[:, 0])**2 + (coords_A_np_arr[:, 1, np.newaxis] - coords_B_np_arr[:, 1])**2 + (coords_A_np_arr[:, 2, np.newaxis] - coords_B_np_arr[:, 2])**2)
    mindist = np.min(dist_bw_atoms_subsystems)
    # minid = np.argmin(dist_bw_atoms_subsystems)
    # minid = np.where(dist_bw_atoms_subsystems == np.min(dist_bw_atoms_subsystems))
    minid = divmod(dist_bw_atoms_subsystems.argmin(), dist_bw_atoms_subsystems.shape[1])
    st.info('Distance b/w the COMs of the subsystems is ' + str(np.round(dist_bw_COM_subsystems, 6))+'  Angstroms', icon='✨')
    st.info('Minimum distance b/w the atoms of the subsystems is ' + str(np.round(mindist, 5))+'  Angstroms between atoms '+ str(selected_rows_A.iloc[minid[0],0])+'('+str(minid[0]+1)+')'+' and '+ str(selected_rows_B.iloc[minid[1],0])+'('+str(minid[1]+1)+')', icon='✨')


    col1, col2 = st.columns(2)
    with col1:
        st.write('#### Selected atoms for subsystem A')
        st.dataframe(selected_rows_A, height=350)
        st.write('Center of Mass of the subsystem A', com_A)
        st.write('No. of atoms in subsystem A', natoms_A)
    with col2:
        st.write('#### Selected atoms for subsystem B')
        st.dataframe(selected_rows_B, height=350)
        st.write('Center of Mass of the subsystem B', com_B)
        st.write('No. of atoms in subsystem B', natoms_B)



    st.write('#### Reformatted Coordinate files with the atoms belonging to the subystem A in the beginning')
    # neater modified xyz
    mol_final = pybel.readstring('xyz', modified_xyz)
    modified_xyz = mol_final.write('xyz')
    modified_coords_file = mol_final.write('tmol')
    col1_coords, col2_coords = st.columns(2)
    col1_coords.text_area(label='XYZ file contents with subsystem A in the beginning', value = modified_xyz, placeholder = 'Put your text here', height=250, key = 'output_xyz_text_area')
    col2_coords.text_area(label='Coords file contents with subsystem A in the beginning', value = modified_coords_file, placeholder = 'Put your text here', height=250, key = 'output_coords_text_area')
    col1_coords.download_button(
        label="Download the XYZ file for your use (This file is not needed for embedding calculations with RIPER)",
        data=modified_xyz,
        file_name='total.xyz',
        mime='text/csv',
    )
    col2_coords.download_button(
        label="Download the totalCoords file for embedding calculations with RIPER",
        data=modified_coords_file,
        file_name='totalCoords',
        mime='text/csv',
    )

    ### BEGIN INPUT FILE CREATION ####
    st.write('#### Set the parameters for the input file')
    ### PARAMETERS ###
    nsystm = 2 # no. of subsystems (For now only two are supported so no option to change this)
    ptnIndx = [natoms_A, natoms_tot] # Indices of the last atoms of the two subssytems separated by space
    KEfunc = 521 # LibXC code of the KEDF (Default: LC94)
    xName = 1 # LibXC code of the exchange functional (Default: Slater Exchange)
    cName = 7 # LibXC code of the correlation functional (Default: VWN5)
    xc_change = False # If a different XC functional should be used for active and encironment subsystems
    xName_cluster = 101 # LibXC code of the exchange functional (Default: PBE)
    cName_cluster = 130 # LibXC code of the exchange functional (Default: PBE)
    basis = 'def2-SVP'  # This is treated as the basis set for the whole system if 'basis_cluster' is not provided 
                        #  in the input file, and the basis set for the environment if the 'basis_cluster' keyword is provided separately
    basis_cluster = 'def2-SVP'
    auxbasis = 'universal'
    #cbas = 'cc-pVDZ'
    method = 1 # Emedding Method
    super = False # If a supermolecular basis is used or not
    charge_cluster = 0
    charge_environment = 0
    frozen = True # If False, a Freeze-and-Thaw (FaT) procedure is performed
    scratch = True # If True, the FaT procedure is performed from scratch
    emb_error = True # If True, a total DFT calculation is performed to calculate the embedding error
    periodicity = 0 # Periodicity for periodic0in-periodic embedding calculations (only relevant for method=5 for now)
    cell_params = '' #Cell parameters for periodic-in-periodic embedding calculations, example: '20.0 20.0 20.0 90.0 90.0 90.0'
    latt_vec_a = ''#Lattice vector a for periodic-in-periodic embedding calculations, example:'20.0 0.0 0.0'
    latt_vec_b = ''#Lattice vector b for periodic-in-periodic embedding calculations, example:'0.0 20.0 0.0'
    latt_vec_c = ''#Lattice vector c for periodic-in-periodic embedding calculations, example:'0.0 0.0 20.0'
    kpoints = '' # k-points grid size for periodic-in-periodic embedding calculations, example: nkpoints 4 4 4
    nmax_FaT = 5 # No. of max FaT cycles if convergence is not reached before
    scfconv = 7 # Energy convergence criteria for both Embedding run and FaT 
    denconv = 7 # Density convergence criteria for both Embedding run and FaT 
    ricore = 500 # Memory to keep RI integrals in core
    mxitdiis = 5 # No. of DIIS vectors to be used
    scfiterlimit = 50 # Max SCF iterations for each Embedding run
    # path = '/home/user/turbomole/bin/sysname/'

    





    ### BASIS SET ####
    st.write('##### Basis related parameters')
    is_same_basis_set = st.checkbox('Use same basis set for the total system', value=True)

    basis_set_list = ('sto-3g', 'def2-SVP', 'def2-SVPD', 'def2-TZVP', 'def2-TZVPP', 'def2-TZVPD', 'def2-TZVPPD', 'def2-QZVP', '6-31G', '6-311G')
    auxbasis_set_list = ('def2-SVP','def2-TZVP','universal')

    if is_same_basis_set:
        basis_set_tot = st.selectbox('Select a basis set for the total system',
            basis_set_list, key='basis_set_tot')
    else:
        basis_set_A = st.selectbox('Select a basis set for the subsystem A',
        basis_set_list, key='basis_set_A')
        basis_set_B = st.selectbox('Select a basis set for the subsystem B',
        basis_set_list, key='basis_set_B')
    
    auxbasis_set =  st.selectbox('Select an auxiliary basis set for the subsystems', auxbasis_set_list, key='auxbasis_set')
    isSuperBasis = st.checkbox('Use a supermolecular basis for the subsystems', value=False)

    ### RI core ####
    st.write('##### RI core')
    ricore_memory = st.number_input('Specify memory size in MB to keep RI integrals in RAM', min_value=50, value=500, step=500)

    ### SCF Convergence related parameters ####
    st.write('##### SCF Convergence related parameters')
    energy_conv = st.number_input('Specify energy convergence criterion (10^-N) for both embedding and Freeze-and-Thaw run ', min_value=3, max_value=10, value=7, step=1)
    density_conv = st.number_input('Specify density convergence criterion (10^-N) for both embedding and Freeze-and-Thaw run ', min_value=3, max_value=10, value=7, step=1)
    max_scf_cycles = st.number_input('Specify the maximum no. of SCF cycles to be performed for each embedding run', min_value=3, max_value=500, value=40, step=1)
    max_it_diis = st.number_input('No. of error vectors for DIIS', min_value=1, max_value=100, value=5, step=1)

    ### Freeze-and-Thaw (FaT) parameters ####
    st.write('##### Freeze-and-Thaw (FaT) parameters')
    isFaT = st.checkbox('Perform Freeze-and-Thaw procedure', value=False)
    if isFaT:
        max_fat_cycles = st.number_input('Specify maximum no. of FaT cycles to be performed if convergence is not reached', min_value=1, max_value=50, value=3, step=1)

    ### Embedding Method ####
    st.write('##### Embedding method')
    embedding_method_list = [1,3,4,5]
    method_code = st.selectbox('Specify the embedding method to be used', embedding_method_list, key='embedding_method')
    method1_description = '''
    ➡️ DFT-in-DFT (molecule-in-molecule) embedding using a KEDF-based embedding potential.\n
    ➡️ Compatible with Freeze-and-Thaw (FaT).\n
    ➡️ One also needs to specify a KEDF for the non-additive Kinetic potential in this case.\n
    ➡️ Gives approximate results (mHa accuracy) for weakly interacting systems even with FaT and supermolecular basis.\n
    ➡️ Large errors for strongly interacting systems.\n
    ➡️ Also, allows to use a supermolecular basis for the subsystems.\n'''
    method3_description = '''
    ➡️ DFT-in-DFT (molecule-in-molecule) embedding using a level-shift projection operator-based embedding potential.\n
    ➡️ Compatible with Freeze-and-Thaw (FaT).\n
    ➡️ Gives exact results, provided a sueprmolecular basis and FaT procedure is used.\n'''
    method4_description = '''
    ➡️ When this method is provided in the control file, then the riper program will read an embedding potential
    from a file 'embdpot' and use it as a fixed additional potential during the SCF iterations. Using the `riperembed.py` script
    or this input file creator is pointless with this method.\n
    ➡️ The expected usage for method 4, is to run an embedding calculation using method=1,2,3,5 and then take the embedding potential saved on disk in a file called 
    'embdpot' and paste it in a directoy where you would like to use this embedding potential as a frozen term for the active subsystem. For example to perform 
    RT-TDDFT calculation on an embedded subsystem using a frozen embedding potential from a molecule-in-molecule(periodic) embedding run.'''
    method5_description = '''
    ➡️ DFT-in-DFT (periodic-in-periodic) embedding using a level-shift projection operator-based embedding potential.\n
    ➡️ Compatible with Freeze-and-Thaw (FaT).\n
    ➡️ Only compatible with a supermolecular basis.\n
    ➡️ Gives exact results.\n
    ➡️ Compatible with multiple k-points.\n
    ➡️ You also need to provide periodicity details like lattice parameters and kpoint grids with this method.\n
    ➡️ Both subsystems are assumed to have the same periodicity, lattice vectors/cell parameters and kmesh.\n'''
    embedding_method_descriptions = {1: method1_description, 3: method3_description, 4: method4_description, 5: method5_description}
    st.write('Description of the chosen embedding method')
    st.write(embedding_method_descriptions[method_code])

    ### Exchange-Correlation Functionals ####
    st.write('##### Exchange-Correlation Functionals')
    is_same_xc = st.checkbox('Use same exchange-correlation functional for both the subsystems', value=True)
    st.write('Use this link to find out more LibXC codes and their references: [https://tddft.org/programs/libxc/functionals/](https://tddft.org/programs/libxc/functionals/)')

    xfunc_dict = {1:'Slater exchange', 101: 'PBE Exchange'}
    cfunc_dict = {7:'VWN5 Correlation', 12: 'Perdew & Wang Correlation', 130: 'PBE Correlation'}
    if is_same_xc:
        xfunc_tot = st.selectbox('Select the LibXC code for exchange functional for the total system',
            xfunc_dict, key='xfunc_tot')
        cfunc_tot = st.selectbox('Select the LibXC code for correlation functional for the total system',
            cfunc_dict, key='cfunc_tot')

    ### Kinetic Energy Functionals ####
    if method_code==1 or method_code==3:
        st.write('##### Kinetic Energy Density Functional')
        st.write('Use this link to find out more LibXC codes and their references: [https://tddft.org/programs/libxc/functionals/](https://tddft.org/programs/libxc/functionals/)')
        kedfunc_dict = {'electro':'none', 521:'LC94 (GGA)', 50: 'Thomas-Fermi KE (LDA)', 55: 'REVAPBE - revised APBE (GGA)', 53: 'REVAPBEINT - interpolated version of revAPBE (GGA)'}
        
        kedfunc = st.selectbox('Select the LibXC code for kinetic energy density functional for the embedding calculations',
                kedfunc_dict, key='kedfunc')

    ### Charges for subsystems ####
    st.write('##### Charges for subsystems')
    st.write('Since only closed-shell occupations are supported for embedding calculations, one may need to provide charges to subsystems to ensure this, especially when covalent bonds are broken.')
    charge_A = st.number_input('Charge for subsystem A', value=0, step=1)
    charge_B = st.number_input('Charge for subsystem B', value=0, step=1)

    ### Periodic details ####
    if method_code==5:
        st.write('##### Details for the periodic system')
        periodicity = st.number_input('Periodicity', value=1, min_value=1, max_value=3, step=1)
        cell_info_mode = st.radio('Will you provide cell parameters (also known as lattice constants/lattice parameters) or lattice vectors?', ('Cell Parameters', 'Lattice Vectors'))
        if cell_info_mode=='Cell Parameters':
            col1_cell, col2_cell, col3_cell = st.columns(3)
            if periodicity==1:
                cell_a = col1_cell.number_input('Cell parameter a (in Angstroms)', value=4.0, step=0.1)
            if periodicity==2:
                cell_a = col1_cell.number_input('Cell parameter a (in Angstroms)', value=4.0, step=0.1)
                cell_b = col2_cell.number_input('Cell parameter b (in Angstroms)', value=4.0, step=0.1)
                cell_gamma = col1_cell.number_input('Cell parameter gamma (in degrees)', value=90.0, step=0.5)
            if periodicity==3:
                cell_a = col1_cell.number_input('Cell parameter a (in Angstroms)', value=4.0, step=0.1)
                cell_b = col2_cell.number_input('Cell parameter b (in Angstroms)', value=4.0, step=0.1)
                cell_c = col3_cell.number_input('Cell parameter c (in Angstroms)', value=4.0, step=0.1)
                cell_alpha = col1_cell.number_input('Cell parameter alpha (in degrees)', value=90.0, step=0.5)
                cell_beta = col2_cell.number_input('Cell parameter beta (in degrees)', value=90.0, step=0.5)
                cell_gamma = col3_cell.number_input('Cell parameter gamma (in degrees)', value=90.0, step=0.5)
        elif cell_info_mode=='Lattice Vectors':
            col1_cell, col2_cell, col3_cell = st.columns(3)
            if periodicity==1:
                x_latt_vec_a = col1_cell.number_input('X-component of Lattice vector a (in Angstroms)', value=4.0, step=0.1)
            if periodicity==2:
                x_latt_vec_a = col1_cell.number_input('X-component of Lattice vector a (in Angstroms)', value=4.0, step=0.1)
                y_latt_vec_a = col2_cell.number_input('Y-component of Lattice vector a (in Angstroms)', value=0.0, step=0.1)
                x_latt_vec_b = col1_cell.number_input('X-component of Lattice vector b (in Angstroms)', value=0.0, step=0.1)
                y_latt_vec_b = col2_cell.number_input('Y-component of Lattice vector b (in Angstroms)', value=4.0, step=0.1)
            if periodicity==3:
                x_latt_vec_a = col1_cell.number_input('X-component of Lattice vector a (in Angstroms)', value=4.0, step=0.1)
                y_latt_vec_a = col2_cell.number_input('Y-component of Lattice vector a (in Angstroms)', value=0.0, step=0.1)
                z_latt_vec_a = col3_cell.number_input('Z-component of Lattice vector a (in Angstroms)', value=0.0, step=0.1)
                x_latt_vec_b = col1_cell.number_input('X-component of Lattice vector b (in Angstroms)', value=0.0, step=0.1)
                y_latt_vec_b = col2_cell.number_input('Y-component of Lattice vector b (in Angstroms)', value=4.0, step=0.1)
                z_latt_vec_b = col3_cell.number_input('Z-component of Lattice vector b (in Angstroms)', value=0.0, step=0.1)
                x_latt_vec_c = col1_cell.number_input('X-component of Lattice vector c (in Angstroms)', value=0.0, step=0.1)
                y_latt_vec_c = col2_cell.number_input('Y-component of Lattice vector c (in Angstroms)', value=0.0, step=0.1)
                z_latt_vec_c = col3_cell.number_input('Z-component of Lattice vector c (in Angstroms)', value=4.0, step=0.1)

        st.write('###### K-mesh/gridsize info')
        col1_kpts, col2_kpts, col3_kpts = st.columns(3)
        if periodicity==1:
            nk_x = col1_kpts.number_input('No. of k-points along x', value=3, step=1)
        if periodicity==2:
            nk_x = col1_kpts.number_input('No. of k-points along x', value=3, step=1)
            nk_y = col2_kpts.number_input('No. of k-points along y', value=3, step=1)
        if periodicity==3:
            nk_x = col1_kpts.number_input('No. of k-points along x', value=3, step=1)
            nk_y = col2_kpts.number_input('No. of k-points along y', value=3, step=1)
            nk_z = col3_kpts.number_input('No. of k-points along z', value=3, step=1)

    ### Calculate Embedding Error ####
    st.write('##### Calculate embedding error?')
    isEmbError = st.checkbox('Should a total regular KS-DFT calculation be performed at the end to estimate the embedding error?', value=True)

    ### RIPER path ####
    # riper_path = st.text_input('Path of the riper or riper_smp/riper_omp executable', value='/home/user/turbomole/bin/em64t-unknown-linux-gnu_smp/')

    ### Start creating the text for the input file ####
    st.write('#### INPUT FILE')
    input_file_str = '# INPUT FILE FOR RUNNING EMBEDDING CALCULATIONS VIA riperembed.py SCRIPT AND RIPER MODULE OF TURBOMOLE\n'
    input_file_str = input_file_str + '''# Cite the implementation as: 
# Manas Sharma and Marek Sierka
# Journal of Chemical Theory and Computation 2022 18 (11), 6892-6904
# DOI: 10.1021/acs.jctc.2c00380\n'''
    input_file_str = input_file_str + '$FDE'
    input_file_str = input_file_str + '\nnsystm = ' + str(nsystm)
    input_file_str = input_file_str + '\nptnIndx = ' + str(ptnIndx[0]) + ' ' + str(ptnIndx[1])
    input_file_str = input_file_str + '\nKEfunc = ' + str(kedfunc)
    input_file_str = input_file_str + '\nbasis = ' + str(basis_set_tot)
    input_file_str = input_file_str + '\nauxbasis = ' + str(auxbasis_set)
    input_file_str = input_file_str + '\nKEfunc = ' + str(kedfunc)
    input_file_str = input_file_str + '\nxName = ' + str(xfunc_tot)
    input_file_str = input_file_str + '\ncName = ' + str(cfunc_tot)
    input_file_str = input_file_str + '\nfrozen = ' + str(not isFaT)
    input_file_str = input_file_str + '\nscratch = ' + str(True)
    input_file_str = input_file_str + '\nsuper = ' + str(isSuperBasis)
    input_file_str = input_file_str + '\nemb_error = ' + str(isEmbError)
    input_file_str = input_file_str + '\nmethod = ' + str(method_code)
    input_file_str = input_file_str + '\nscfconv = ' + str(energy_conv)
    input_file_str = input_file_str + '\ndenconv = ' + str(density_conv)
    input_file_str = input_file_str + '\nricore = ' + str(ricore_memory)
    if isFaT:
        input_file_str = input_file_str + '\nnmax_FaT = ' + str(max_fat_cycles)
    input_file_str = input_file_str + '\nscfiterlimit = ' + str(max_scf_cycles)
    input_file_str = input_file_str + '\nmxitdiis = ' + str(max_it_diis)
    if method_code==5:
        input_file_str = input_file_str + '\nperiodicity = ' + str(periodicity)
    # input_file_str = input_file_str + '\npath = ' + str(riper_path)



    st.text_area('Input file contents for embedding via RIPER', value=input_file_str, height=400)
    st.download_button(
        label="Download the input file",
        data=input_file_str,
        file_name='input',
        mime='text/csv',
    )

    st.write('#### What next?')
    st.write('1️⃣ Download the `totalCoords` file from the previous section.')
    st.write('2️⃣ Download the `input` file that we just created.')
    st.write('3️⃣ Put the two files in the same directory.')
    st.write('4️⃣ Run the `riperembed.py` script as: `nohup riperembed.py > output_embedding &`')

    st.write('### Further information to understand the outputs')
    st.write(further_information)
