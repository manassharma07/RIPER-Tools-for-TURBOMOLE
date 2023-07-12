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



st.title("Relevant keywords for RIPER")


st.write("""`riper` shares most of the relevant keywords of the `dscf` and `ridft` modules. The
`$dft` data group and auxiliary basis sets defined using the keyword
`$jbas` are always required.
For periodic calculations two additional keywords are necessary:
""")
data_periodic_keyword = {
    'Keyword': [
        '$periodic n',
        '$cell', 
        '$lattice',
    ],
    'Meaning': [
        'Specifies the number of periodic directions: n = 3 for a 3D periodic bulk solid, n = 2 for a 2D periodic surface slab and n = 1 for a 1D periodic system. The default value is 0 for a molecular system.',
        'Specifies the unit cell parameters. The number of cell parameters depends on the periodicity of the system: For 3D periodic systems six unit cell parameters |a|, |b|, |c|, α, β and γ need to be provided. Here, |a|, |b| and |c| are lengths of the appropriate cell vectors, α is the angle between vectors b and c, β is the angle between vectors a and c, and γ is the angle between vectors a and b. riper assumes that the cell vectors a and b are aligned along the x axis and on the xy plane, respectively. For 2D periodic systems three surface cell parameters |a|, |b| and γ have to be provided. Here, |a| and |b| are lengths of the appropriate cell vectors and γ is the angle between a and b. riper assumes that the cell vectors a and b are aligned along the x axis and on the xy plane, respectively. For 1D periodic systems only one parameter specifying the length of the unit cell has to be provided. riper assumes that periodic direction is along the x axis.',
        '''Alternatively, lattice vectors can be provided. The number of cell parameters
depends on the periodicity of the system:
For 3D periodic systems three (three-dimensional) lattice vectors need to be
provided.
For 2D periodic systems two (two-dimensional) lattice vectors have to be provided. riper assumes that the lattice vectors are aligned on the xy plane.
For 1D periodic systems only one parameter specifying the length of the lattice
vector has to be provided. riper assumes that periodic direction is along the
x axis.
        ''',
    ]
}
df_periodic_keywords = pd.DataFrame(data_periodic_keyword)
st.table(df_periodic_keywords)

with st.expander(label='#### `$riper` keywords', expanded=True):
    st.write("""Keywords within the section `$riper` can be used to control precision and 
            parameters of algorithms implemented in riper. If the `$riper` group is absent, 
            the following default values are used""")

    code_riper_keywords = """
    $riper
        # general keywords
        mxitdiis 5
        thrints 1.0d-12
        lenonly off
        lchgprj on (for periodic systems)
        lchgprj off (for molecular systems)
        northol 5
        pqmatdiag off
        pqsingtol 1.0d-8
        # CFMM control options
        lmaxmom 20
        nctrgt 10 (for periodic systems)
        nctrgt 1 (for molecular systems)
        wsicl 3.0
        epsbext 1.0d-9
        locmult on (for periodic systems)
        locmult off (for molecular systems)
        locmomadd 2
        # LMIDF control options
        lpcg on
        lcfmmpcg on
        lmxmompcg 20
        pcgtol 1.0d-9
        pcgtyp sp
        # Gaussian smearing options
        sigma 0.0d0
        desnue 0.0d0
    """
    st.code(code_riper_keywords, language="shell")

    data_riper_keywords = {
        'Keyword': [
            'mxitdiis',
            'thrints',
            'lenonly',
            'lchgprj',
            'northol',
            'pqmatdiag',
            'pqsingtol',
            'lmaxmom',
            'nctrgt',
            'wsicl',
            'epsbext',
            'locmult',
            'locmomadd',
            'lpcg',
            'lcfmmpcg',
            'lmxmompcg',
            'pcgtol',
            'pcgtyp',
            'sigma',
            'desnue'
        ],
        'Meaning': [
            'DIIS subspace size (number of previous KS/Fock matrices used for determining the DIIS coefficients)',
            'Threshold for integrals neglect and for differential overlap when screening basis functions pairs. Probably never needs to be changed.',
            'Flag for energy calculation only, no gradients.',
            'If set to on charge projection of the auxiliary electron density [164] is performed for molecular systems during calculation of the Coulomb term. The charge projection constraints the charge of auxiliary density exactly to the number of electrons in the system. It is required for periodic systems, otherwise the Coulomb energy would be infinitely large. For molecular systems charge projection leads to a slight increase of the RI fitting error. It may be useful in some cases but we have so far not identified any.',
            'Forces orthonormalization of orbital coefficients every northol SCF iteration.',
            'If set to on full diagonalization of the Coulomb metric matrix [164] is performed and used to solve density fitting equations. When diffuse auxiliary basis functions are used the default Cholesky decomposition of the Coulomb metric matrix may fail due to small negative eigenvalues. In this case the slower method based on a full diagonalization of the metric matrix is necessary.',
            'If pqmatdiag is used pqsingtol sets threshold for neglect of small eigenvalues of the Coulomb metric matrix.',
            'Maximum l-moment of multipole expansions used for calculation of the Coulomb term. The default value hardly ever needs to be changed.',
            'Target number of charge distributions per lowest level box of the octree [160]. The default value hardly ever needs to be changed.',
            'Sets the well-separateness criterion [160]. Octree boxes with centers separated more than sum of their lengths times 0.5×wsicl are considered as well-separated. The default hardly ever needs to be changed. wsicl makes sense only for values ≥ 2.0. For wsicl ≤ 3.0 increasing lmaxmom may be necessary for reasonable accuracy.',
            'Precision parameter used to determine basis function extents [160].',
            'If set to on, an additional acceleration method employing local multipole expansions is used. For periodic systems this leads to a significant speedup of calculations, especially for small unit cells and/or diffuse basis functions. Default value is off and on for molecular and periodic systems, respectively.',
            'For locmult set to on the order of local multipole expansions is increased by locmomadd. The default value probably never needs to be changed.',
            'If set to on the low-memory iterative density fitting (LMIDF) scheme is used for solving the RI equations [163] using the preconditioned conjugate gradient (PCG) method. It is implemented for molecular systems only. Default value is off.',
            'If lpcg is used, lcfmmpcg specifies whether the CFMM is applied for evaluation of the matrix-vector products needed for the PCG solver. Not employing CFMM speeds up the calculations but significantly increases memory demand since the full Coulomb metric matrix has to be stored. Default value is on.',
            'Maximum l-moment of multipole expansions for calculations of Coulomb interactions within the PCG algorithm. It should be set to the same or larger value than lmaxmom.',
            'Sets the threshold parameter controlling accuracy of the PCG solver (see [163] for details). Default value is 1.0 · 10^−9. For lower-precision calculations it can be set to 1.0 · 10^−8 but values larger than 1.0 · 10^−7 are not allowed as these lead to large errors in Coulomb energies and occasionally to SCF convergence problems.',
            'Specifies the type of preconditioner used in the PCG algorithm. Three types of preconditioners are implemented and are defined explicitly in Sec. 7. The sp preconditioner is a default one performing consistently the best among the preconditioners considered. The at preconditioner is less efficient in decreasing the number of CG iterations needed for convergence. However, it has negligible memory requirements and more favorable scaling properties, albeit with a large prefactor. The ss preconditioner represents a middle ground between the sp and at preconditioners both in terms of the efficiency and memory requirements.',
            'Width of the Gaussian smearing in hartree. See ref. [166] for more information. Note that [166] uses eV as the unit for the width of the Gaussian smearing.',
            'Specifying desnue along with sigma forces occupancy leading to the number of unpaired electrons equal to desnue.'
        ]
    }

    df_riper_keywords = pd.DataFrame(data_riper_keywords)
    st.table(df_riper_keywords)