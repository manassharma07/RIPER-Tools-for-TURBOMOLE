import streamlit as st
import webbrowser
webbrowser.open_new_tab('https://cubesuite.streamlit.app')

# Set page config
st.set_page_config(page_title='Cube File Processing (Gateway)', layout='wide', page_icon="⚛️",
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


st.write('# Cube File Processing')
st.write('###### You can process cube files used to plot electronic density & orbitals by visiting the ✨[**Cube Suite**](https://cubesuite.streamlit.app)✨ web app from the [**Riper-Tools**](https://ripertools.turbomole.org) developer.')

st.page_link('https://cubesuite.streamlit.app', label='# 🔗 Visit **Cube Suite**')

st.write('###### ➡️[**Cube Suite**](https://cubesuite.streamlit.app) is powered by the amazing [Cube Toolz](https://github.com/funkymunkycool/Cube-Toolz) python package under the hood and is basically a GUI for [Cube Toolz](https://github.com/funkymunkycool/Cube-Toolz).')
st.write('###### [**Cube Suite**](https://cubesuite.streamlit.app) can help you perform tasks like:')
st.write('    ➕ [Add two Cube files](https://cubesuite.streamlit.app/Add_Cube_Files)')
st.write('    ➖ [Subtract two Cube files](https://cubesuite.streamlit.app/Subtract_Cube_Files)')
st.write('    ✖️ [Multiply two Cube files](https://cubesuite.streamlit.app/Multiply_Cube_Files)')
st.write('    ↗️ [Translate a Cube file](https://cubesuite.streamlit.app/Translate_Cube_File)')
st.write('    📉 [Calculate Planar Average of a Cube file](https://cubesuite.streamlit.app/Planar_Average_of_Cube_File)')


