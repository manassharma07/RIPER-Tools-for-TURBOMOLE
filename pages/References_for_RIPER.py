import streamlit as st
import pubchempy as pcp
from pymatgen.core.structure import Molecule
import py3Dmol
import pandas as pd
import streamlit.components.v1 as components

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
st.sidebar.write('### *Contributors*')
st.sidebar.write('[Ya-Fan Chen ](https://github.com/Lexachoc)')
st.sidebar.write('### *Source Code*')
st.sidebar.write('[GitHub Repository](https://github.com/manassharma07/RIPER-Tools-for-TURBOMOLE)')



st.title("Relevant references for RIPER")



references = [
    {
        'Author': 'Burow, A.M.; Sierka, M.; Mohamed, F.',
        'Title': 'Resolution of identity approximation for the Coulomb term in molecular and periodic systems.',
        'Journal': 'J. Chem. Phys. 2009, 131, 214101-1-214101-6. (56)',
        'URL': 'https://doi.org/10.1063/1.3267858'
    },
    {
        'Author': 'Łazarski, R.; Burow, A.M.; Sierka, M.',
        'Title': 'Density functional theory for molecular and periodic systems using density fitting and continuous fast multipole methods.',
        'Journal': 'J. Chem. Theory Comput. 2015, 11, 3029-3041. (57)',
        'URL': 'https://pubs.acs.org/doi/10.1021/acs.jctc.5b00252'
    },
    {
        'Author': 'Łazarski, R.; Burow, A.M.; Grajciar, L.; Sierka, M.',
        'Title': 'Density functional theory for molecular and periodic systems using density fitting and continuous fast multipole method: Analytical gradients.',
        'Journal': 'J. Comput. Chem. 2016, 37, 2518-2526. (58)',
        'URL': 'https://doi.org/10.1002/jcc.24477'
    },
    {
        'Author': 'Becker, M.; Sierka, M.',
        'Title': 'Density functional theory for molecular and periodic systems using density fitting and continuous fast multipole method: Stress tensor.',
        'Journal': 'J. Comput. Chem. 2019, 40, 2563-2570.',
        'URL': 'https://doi.org/10.1002/jcc.26033'
    },
    {
        'Author': 'Grajciar, L.',
        'Title': 'Low-memory Iterative Density Fitting.',
        'Journal': 'J. Comput. Chem. 2015, 36, 1521-1535. (327)',
        'URL': 'https://doi.org/10.1002/jcc.23961'
    },
    {
        'Author': 'Burow, A.M.; Sierka, M.',
        'Title': 'Linear scaling hierarchical integration scheme for the exchange-correlation term in molecular and periodic systems.',
        'Journal': 'J. Chem. Theory Comput. 2011, 7, 3097-3104. (328)',
        'URL': 'https://pubs.acs.org/doi/10.1021/ct200412r'
    },
    {
        'Author': 'Irmler, A.; Burow, A.M.; Pauly, F.',
        'Title': 'Robust Periodic Fock Exchange with Atom-Centered Gaussian Basis Sets.',
        'Journal': 'J. Chem. Theory Comput. 2018, 14, 4567-4580.',
        'URL': 'https://pubs.acs.org/doi/full/10.1021/acs.jctc.8b00122'
    },
    {
        'Author': 'Müller, C.; Sharma, M.; Sierka, M.',
        'Title': 'Real-time time-dependent density functional theory using density fitting and the continuous fast multipole method.',
        'Journal': 'J. Comput. Chem. 2020, 41, 2573-2582.',
        'URL': 'https://doi.org/10.1002/jcc.26412'
    },
    {
        'Author': 'Sharma, M.; Sierka, M.',
        'Title': 'Efficient Implementation of Density Functional Theory Based Embedding for Molecular and Periodic Systems Using Gaussian Basis Functions.',
        'Journal': 'J. Chem. Theory Comput. 2022, 18, 6892-6904',
        'URL': 'https://pubs.acs.org/doi/10.1021/acs.jctc.2c00380'
    }
]

df = pd.DataFrame(references)

# Format the 'Author' column
df['Author'] = df['Author'].apply(lambda author: author.replace(';', ', '))

# Format the 'Journal' column with hyperlinks
df['Journal'] = df.apply(lambda row: f'<a href="{row["URL"]}">{row["Journal"]}</a>', axis=1)

# Drop the 'URL' column
df.drop('URL', axis=1, inplace=True)

# Display the references as a table with clickable links
df = df.to_html(escape=False)
st.write(df, unsafe_allow_html=True)