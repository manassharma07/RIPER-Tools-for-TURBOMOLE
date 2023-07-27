import streamlit as st
import pubchempy as pcp
from pymatgen.core.structure import Molecule
import py3Dmol
import pandas as pd
import streamlit.components.v1 as components

# Set page config
st.set_page_config(page_title='Convergence Tips for RIPER', layout='wide', page_icon="⚛️",
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



st.title("Convergence tips for RIPER")


st.write('When using `riper` for calculations on metals, and facing convergence issues, one can try the following settings:')
st.write('Increase the number of DIIS matrices and the accuracy of lattice sums')
code_riper_metals = """
    $riper
        sigma 0.01
        epsbext 1.0d-9
        mxitdiis 100
    """
st.code(code_riper_metals, language="shell")
st.write('Increase automatic shift and damping')
code_riper_metals = """
    $scforbitalshift  automatic=1.0
    $scfdamp   start=5.00  step=0.050  min=0.50
    """
st.code(code_riper_metals, language="shell")

  