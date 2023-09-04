import streamlit as st
import platform

# Set page config
st.set_page_config(page_title='RIPER Tools for TURBOMOLE', layout='wide', page_icon="⚛️",
menu_items={
         'About': "A web app to help you with DFT related calculations using the RIPER module of [TURBOMOLE](https://www.turbomole.org/)"
     })

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


# Main app
st.header('RIPER Tools')
st.write('#### A web app to help you with DFT related calculations using the RIPER module of [TURBOMOLE](https://www.turbomole.org/)')

st.write('''The riper module is an implementation of Kohn-Sham DFT with Gaussian-type orbitals (GTO) as basis functions that treats molecular and periodic systems of any dimensionality on an equal footing.''')
# "to create input (control) files from a given [CIF](https://en.wikipedia.org/wiki/Crystallographic_Information_File), [XYZ](https://en.wikipedia.org/wiki/XYZ_file_format), etc"
# "You can also parse the output of the RIPER program for further analysis."
text_intro = """
This app allows you to create `RIPER` input files for a material in the [materialsproject.org](https://next-gen.materialsproject.org/) database or a molecule in the [PubChem database](https://pubchem.ncbi.nlm.nih.gov/). 
Additionally it allows to convert `CIF`s, `XYZ`s, `POSCAR`s and other popular formats to `RIPER` format.
You can also use it to get help with `RIPER` related examples, keywords, references, and tutorials. 
The web app can also do things like band structure path generation for band structure calculations, plot density of states, generate input for RT-TDDFT simulations, parse output of `RIPER` calculations, and perform unit conversions.

You can select from a variety of options from the menu in the left sidebar."""
st.write(text_intro)

st.write('By default `riper` calculates single point energy and gradient in one run. For this, simply invoke riper as')
st.code('nohup riper > riper.out &', language='shell')

st.write('''For energy calculation only add `lenonly` on to the `$riper` section in the `control` file, i.e.''')
st.code('''
$riper
    lenonly on
''', language='shell')

st.write('#### Tutorial #1 (Molecular DFT calculation using RIPER)')
_, container, _ = st.columns([50, 100, 50])
container.video(data='https://www.youtube.com/watch?v=5ZlbBYxjwX8')

st.write(''' `jobex` will automatically use riper when the keyword `$periodic` is present in the `control`
file. Alternatively, the use of riper can be forced by specifying `-riper` argument of `jobex`,
i.e., invoking''')
st.code('nohup jobex -riper > jobex.out &', language='shell')

st.write('''Simultaneous optimization of atomic positions and cell parameters or lattice vectors can be
performed by specifying the keyword `$optcell` in the `control` file and then invoking the
jobex script as described above.''')

st.write('#### Tutorial #2 (Molecular Geometry Optimization using RIPER)')
_, container, _ = st.columns([50, 100, 50])
container.video(data='https://www.youtube.com/watch?v=-Io3hpYVs84')

st.write('''
Calculation of data on grids for molecular and periodic systems using `riper` requires the
keyword `$pointvalper` in the `control` file. The values are calculated on a 3D grid and
written to appropriate output files. These files can be generated either running a single
point energy riper calculation or invoking
''')
st.code('nohup riper -proper > riper.out &', language='shell')
st.write('#### Tutorial #3 (Calculating Densities and MOs on grids using RIPER)')
_, container, _ = st.columns([50, 100, 50])
container.video(data='https://www.youtube.com/watch?v=549bXQvaAjU')
