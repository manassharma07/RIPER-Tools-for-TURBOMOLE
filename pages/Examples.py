import streamlit as st
import pubchempy as pcp
from pymatgen.core.structure import Molecule
import py3Dmol
import pandas as pd
import streamlit.components.v1 as components

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



st.title("Examples for calculations using RIPER")


st.write('### Example 1')
st.write('1D periodic system with a unit cell |a| = 5 bohr:')
st.code('''
$periodic 1
$cell
    5.0
''', language="shell")
st.write('The same input using the `$lattice` keyword')
st.code('''
$periodic 1
$lattice
    5.000000000000000
''', language="shell")

st.write('### Example 2')
st.write('''2D periodic system with a unit cell |a| = 5 and |b| = 8 bohr, the angle
between a and b of 60 degrees:''')
st.code('''
$periodic 2
$cell
    5.0 8.0 60.0
''', language="shell")

st.write('The same input using the `$lattice` keyword')
st.code('''
$periodic 2
$lattice
    5.000000000000000 0.000000000000000
    6.928203230275509 4.000000000000000
''', language="shell")

st.write('### Example 3')
st.write('''
3D periodic system with a triclinic unit cell |a| = 5, |b| = 8 and |c| = 6
bohr, the angle between b and c of 70 degrees, between a and c of 60 degrees, and between a and
b of 90 degrees:
''')
st.code('''
$periodic 3
$cell
    5.0 8.0 6.0 70.0 60.0 90.0

''', language="shell")
st.write('The same input using the `$lattice` keyword')
st.code('''
$periodic 3
$lattice
    5.000000000000000 0.000000000000000 0.000000000000000
    6.928203230275509 4.000000000000000 0.000000000000000
    4.328764234233443 3.324345345565466 6.354331231232453

''', language="shell")

st.write('### Example 4')
st.write('''
 In this example a 3D periodic system is defined (`$periodic 3`). The unit
cell is specified using the `$cell` keyword, with lengths and angles given in atomic units and
degrees, respectively. A uniform grid of 27 (3 x 3 x 3) k points for numerical integration over
the BZ is specified using the keyword `$kpoints`. The section `$riper` is
added with the keyword `lenonly` set to `on`. It forces riper to skip the calculation of energy
gradients (by default both energy and gradients are calculated in one run). In addition,
Gaussian smearing is switched on by providing `sigma 0.01`.
''')
st.code('''
$periodic 3
$cell
    18.5911 16.5747 16.5747 90. 90. 90.
$kpoints
    nkpoints 3 3 3
$riper
    lenonly on
    sigma 0.01

''', language="shell")

st.write('The same input using the `$lattice` keyword')
st.code('''
$periodic 3
$lattice
    18.5911 0.0000 0.0000
    0.0000 16.5747 0.0000
    0.0000 0.0000 16.5747
$kpoints
    nkpoints 3 3 3
$riper
    lenonly on
    sigma 0.01

''', language="shell")

st.write('### Example 5')
st.write('''
In this example a 2D periodic system is defined (`$periodic 2`). The unit
cell is specified using the `$cell` keyword, with lengths and angles given in Angs and degrees,
respectively. A uniform grid of 9 (3 x 3) k-points for numerical integration over the BZ is
defined using the keyword `$kpoints`. The section `$riper` is added with
the keyword `lmaxmom 30`. This changes the maximum order of multipole expansions used
in CFMM to 30 (default value is 20).
''')
st.code('''
$periodic 2
$cell angs
    18.5911 16.5747 90.0
$kpoints
    nkpoints 3 3
$riper
    lmaxmom 30
''', language="shell")

st.write('The same input using the `$lattice` keyword')
st.code('''
$periodic 2
$lattice angs
    18.5911 0.0000
    0.0000 16.5747
$kpoints
    nkpoints 3 3
$riper
    lmaxmom 30
''', language="shell")


st.write('### Example 6')
st.write('''
In this example a 1D periodic system is defined (`$periodic 1`). The dimension of the periodic direction in atomic units is specified using the `$cell` keyword. 
A one dimensional grid of 3 k-points for numerical integration over the BZ is defined using the
keyword `$kpoints`.

''')
st.code('''
$periodic 1
$cell
    18.5911
$kpoints
    nkpoints 3
''', language="shell")

st.write('The same input using the `$lattice` keyword')
st.code('''
$periodic 1
$lattice
    18.5911
$kpoints
    nkpoints 3

''', language="shell")

st.write('### Example 7')
st.write('''
In this example a molecular system is defined (default if no `$periodic` is
specified). The section `$riper` is added with the keyword `lpcg on`, which
activates the low-memory modification of the RI approximation for calculations on very
large molecular systems
''')
st.code('''
$riper
    lpcg on
''', language="shell")