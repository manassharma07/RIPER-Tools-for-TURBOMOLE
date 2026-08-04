import gzip
import io
import re

import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


# Set page config
st.set_page_config(page_title='Band Structure Plotter for RIPER', layout='wide', page_icon="⚛️",
menu_items={
         'About': "A web app to help you with DFT related calculations using the RIPER module of [TURBOMOLE](https://www.turbomole.org/)"
     })

# Sidebar stuff
st.sidebar.write('# About')
st.sidebar.write(' Originally Made By [Manas Sharma](https://manas.bragitoff.com)')
st.sidebar.write(' In the group of [Prof. Dr. Marek Sierka](https://cmsg.uni-jena.de)')
st.sidebar.write('## Cite us:')
st.sidebar.write('[J. Phys. Chem. A 2025, 129, 39, 9062–9083](https://doi.org/10.1021/acs.jpca.5c02937)')
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


HARTREE2EV = 27.211324570273
PALETTE = ['#e41a1c', '#377eb8', '#4daf4a', '#ff7f00', '#984ea3', '#a65628', '#f781bf',
           '#17becf', '#bcbd22', '#8c564b', '#e377c2', '#7f7f7f']

LOOKS_LIKE = [
    ('$coord', 'this looks like a **coord** file (the structure), not a band structure file'),
    ('$cell', 'this looks like a **coord** file (the structure), not a band structure file'),
    ('$kpoints', 'this looks like the **control** file. It *requests* the band structure; the '
                 'bands themselves are written to `bands.xyz` when you run `riper`'),
    ('$scfmo', 'this looks like an **mos**/orbital file, not a band structure file'),
    ('#     DOS', 'this is a **dos** file. Use the *DOS Plot* page for it'),
    ('#     PDOS', 'this is a **pdos** file. Use the *DOS Plot* page for it'),
]


# ----------------------------------------------------------------------------------------
# Parsing (cached). Heavy lifting happens once per uploaded file: everything downstream
# (widget changes, replotting) reuses the cached arrays, which is what keeps the page fast
# on Streamlit Cloud even for bands.xyz files with hundreds of core bands. The numeric block
# is parsed by pandas' C engine, so no numba/JIT is needed (and cold-start JIT would actually
# be slower on the Cloud).
# ----------------------------------------------------------------------------------------
@st.cache_data(show_spinner='Parsing band structure file ...')
def parse_bands(raw: bytes):
    """Parse bands.xyz / bands_proj.xyz.

    Returns dict with:
      kxyz (npts,3), kpath (npts,), E {spin: (nbands,npts)}, W {spin: (nbands,npts,nw)},
      buckets [labels], npts, nbands, nspin, projected
    Raises ValueError with a user-facing message on unusable input.
    """
    if not raw:
        raise ValueError('the file is empty')
    if raw[:2] == b'\x1f\x8b':          # gzip magic: accept bands.xyz.gz transparently
        try:
            raw = gzip.decompress(raw)
        except Exception:
            raise ValueError('the .gz file could not be decompressed')
    if not raw.strip():
        raise ValueError('the file is empty')
    text = raw.decode('utf-8', errors='replace')

    head = text[:4000]
    for token, hint in LOOKS_LIKE:
        if token in head:
            raise ValueError(hint)

    # bucket labels from the bands_proj.xyz header
    buckets = []
    for m in re.finditer(r'^#\s*bucket\s+(\d+)\s*:\s*(.+)$', head, re.M):
        buckets.append(m.group(2).strip())

    # split off the spin sections; the section markers are the only non-numeric,
    # non-comment lines
    alpha_tag, beta_tag = 'closed/alpha shells:', 'beta shells:'
    ia = text.find(alpha_tag)
    if ia < 0:
        raise ValueError('no `closed/alpha shells:` marker found. Expected a `bands.xyz` or '
                         '`bands_proj.xyz` file as written by riper')
    ib = text.find(beta_tag)
    sections = {}
    if ib < 0:
        sections[1] = text[ia + len(alpha_tag):]
    else:
        sections[1] = text[ia + len(alpha_tag):ib]
        sections[2] = text[ib + len(beta_tag):]

    data = {}
    ncols = None
    for spin, blk in sections.items():
        try:
            df = pd.read_csv(io.StringIO(blk), sep=r'\s+', header=None, comment='#',
                             engine='c', dtype=np.float64)
        except Exception:
            raise ValueError('the numeric rows could not be parsed. The file may be corrupted '
                             'or truncated (e.g. copied while riper was still writing it)')
        if df.shape[1] < 5:
            raise ValueError('rows have %d columns; expected at least 5 '
                             '(kx ky kz |k| energy)' % df.shape[1])
        if ncols is None:
            ncols = df.shape[1]
        elif df.shape[1] != ncols:
            raise ValueError('the alpha and beta sections have different numbers of columns')
        data[spin] = df.to_numpy()

    nw = ncols - 5
    if nw > 0 and len(buckets) != nw:
        buckets = ['bucket %d' % (i + 1) for i in range(nw)]

    arr = data[1]
    # number of k points per band: the k path restarts at the first point for every band
    same = np.flatnonzero((arr[:, 3] == arr[0, 3]) & (arr[:, 0] == arr[0, 0]) &
                          (arr[:, 1] == arr[0, 1]) & (arr[:, 2] == arr[0, 2]))
    npts = int(same[1]) if len(same) > 1 else arr.shape[0]
    if arr.shape[0] % npts != 0:
        raise ValueError('could not infer the number of k points per band (%d rows, first '
                         'repeat at %d). Is the file complete?' % (arr.shape[0], npts))
    nbands = arr.shape[0] // npts

    E, W = {}, {}
    for spin, a in data.items():
        if a.shape[0] != nbands * npts:
            raise ValueError('the beta section has a different size than the alpha section')
        E[spin] = np.ascontiguousarray(a[:, 4].reshape(nbands, npts))
        if nw > 0:
            W[spin] = np.ascontiguousarray(a[:, 5:].reshape(nbands, npts, nw))

    return {'kxyz': data[1][:npts, 0:3].copy(), 'kpath': data[1][:npts, 3].copy(),
            'E': E, 'W': W, 'buckets': buckets, 'npts': npts, 'nbands': nbands,
            'nspin': len(data), 'projected': nw > 0}


@st.cache_data
def detect_ticks(kxyz: np.ndarray, kpath: np.ndarray):
    """Segment boundaries of the k path = direction changes of the k vector."""
    d = np.diff(kxyz, axis=0)
    n = np.linalg.norm(d, axis=1)
    n[n < 1e-30] = 1e-30
    cosang = np.sum(d[1:] * d[:-1], axis=1) / (n[1:] * n[:-1])
    idx = [0] + [j + 1 for j in np.flatnonzero(cosang < 1.0 - 1e-8)] + [len(kpath) - 1]
    idx = sorted(set(idx))
    return [float(kpath[j]) for j in idx]


@st.cache_data
def band_edges(E1: np.ndarray):
    """Per-band extrema and the gap list used for the automatic VBM detection."""
    bmax = E1.max(axis=1)
    bmin = E1.min(axis=1)
    gaps = bmin[1:] - bmax[:-1]          # indirect gap above band i (may be negative)
    return bmax, bmin, gaps


@st.cache_data(show_spinner=False)
def parse_output(raw: bytes):
    """Extract the final Fermi level statistics (a.u.) from a riper output file."""
    if raw[:2] == b'\x1f\x8b':
        try:
            raw = gzip.decompress(raw)
        except Exception:
            raise ValueError('the .gz file could not be decompressed')
    text = raw.decode('utf-8', errors='replace')
    if 'R I P E R' not in text and 'riper' not in text:
        raise ValueError('this does not look like a riper output file')
    blocks = text.split('Fermi Level Statistics (au)')
    if len(blocks) < 2:
        raise ValueError('no `Fermi Level Statistics (au)` block found - did the calculation '
                         'reach the final Fermi statistics?')
    blk = blocks[-1][:1500]
    out = {}
    for key, pat in [('cbm', r'Lowest unoccupied band\s*=\s*(-?\d+\.\d+)'),
                     ('vbm', r'Highest occupied band\s*=\s*(-?\d+\.\d+)'),
                     ('gap', r'Band gap\s*=\s*(-?\d+\.\d+)'),
                     ('fermi', r'Fermi level\s*=\s*(-?\d+\.\d+)')]:
        m = re.search(pat, blk)
        if m:
            out[key] = float(m.group(1))
    if 'fermi' not in out:
        raise ValueError('the Fermi level could not be read from the '
                         '`Fermi Level Statistics (au)` block')
    return out


def gap_candidates(bmax, bmin, gaps):
    """Candidate VBM band indices (0 based) from the inter-band gaps.

    A gap qualifies if it is at least 0.1 eV wide and the band below it disperses by at
    least 0.45 eV over the path. The dispersion test is what separates the fundamental gap
    from the (often much wider) gaps above the flat core and semicore bands - deep cores
    spread by ~0.01 eV, semicore s bands by ~0.3-0.4 eV, while a genuine valence band top
    disperses by 0.5 eV or (usually much) more. The *first* candidate is the best guess:
    holes in the sparse virtual spectrum of a small basis lie above the fundamental gap.
    Returns a list of (index, gap_eV, spread_eV), ordered by band index.
    """
    spread = (bmax - bmin) * HARTREE2EV
    cand = [(i, gaps[i] * HARTREE2EV, spread[i]) for i in range(len(gaps))
            if gaps[i] * HARTREE2EV >= 0.1 and spread[i] >= 0.45]
    return cand


# ----------------------------------------------------------------------------------------
# Page
# ----------------------------------------------------------------------------------------
st.title('Band Structure Plotter')
st.write('Plots `bands.xyz` (plain band structure) and `bands_proj.xyz` (projected / fat '
         'bands from `$pdosper bands on`) written by `riper`. The file type is recognised '
         'automatically. By default the plot zooms in on the bands around the gap instead of '
         'showing hundreds of core bands.')

with st.expander('How to produce these files with RIPER', expanded=False):
    st.write('Band lines are requested inside `$kpoints`; the projection additionally needs '
             '`$pdosper` with `bands on`:')
    st.code("""$kpoints
  nkpoints 4 4 4
  kptlines 3
  recipr 0.0 0.0 0.0  0.5 0.0 0.0  20
  recipr 0.5 0.0 0.0  0.5 0.5 0.0  20
  recipr 0.5 0.5 0.0  0.0 0.0 0.0  21
$pdosper
  elements mg,o     # optional selection: atoms 1,3-5 / resolve atom|l|shell
  bands on""", language='text')
    st.write('`riper` then writes `bands.xyz` and, with `bands on`, `bands_proj.xyz` at the '
             'end of the SCF. Upload either file below.')

uploaded_file = st.file_uploader('Upload `bands.xyz` or `bands_proj.xyz` (gzipped `.gz` also works)')

if uploaded_file is None:
    st.info('⬆️ Upload a `bands.xyz` or `bands_proj.xyz` file to get started. For very large '
            'files, compress first (`gzip bands_proj.xyz`) and upload the `.gz` — band files '
            'typically shrink 5-8x.')
    st.stop()

try:
    B = parse_bands(uploaded_file.getvalue())
except ValueError as exc:
    st.error('**%s** could not be used: %s.' % (uploaded_file.name, str(exc)))
    st.stop()

kind = 'projected band structure (fat bands)' if B['projected'] else 'plain band structure'
st.success('Recognised a **%s**: %d bands × %d k-points%s%s.'
           % (kind, B['nbands'], B['npts'],
              ', spin-polarised (alpha + beta)' if B['nspin'] == 2 else '',
              ', %d weight buckets' % len(B['buckets']) if B['projected'] else ''))

# ----------------------------------------------------------------------------------------
# 1. Fermi level / energy zero
# ----------------------------------------------------------------------------------------
st.subheader('1. Energy reference (Fermi level)')

bmax, bmin, gaps = band_edges(B['E'][1])
cands = gap_candidates(bmax, bmin, gaps)

out_up = st.file_uploader(
    'Optional but recommended: upload the **riper output** file of the same run - it '
    'contains the exact Fermi level, VBM and CBM, so no detection is needed.',
    key='outfile')
outinfo = None
if out_up is not None:
    try:
        outinfo = parse_output(out_up.getvalue())
    except ValueError as exc:
        st.error('**%s** could not be used: %s.' % (out_up.name, str(exc)))

modes = []
if outinfo is not None:
    modes.append('From the riper output (exact)')
modes += ['Automatic: top of the valence band (VBM)',
          'Custom value (e.g. Fermi level from the riper output)']

ef_col1, ef_col2 = st.columns([1.2, 2])
ef_mode = ef_col1.radio('Where to put E = 0?', modes)

if ef_mode.startswith('From the riper output'):
    zero_at = ef_col1.radio('Zero at', ['VBM (highest occupied band)',
                                        'Fermi level (band gap middle)'],
                            disabled='vbm' not in outinfo)
    if zero_at.startswith('VBM') and 'vbm' in outinfo:
        efermi = outinfo['vbm']
    else:
        efermi = outinfo['fermi']
    msg = 'From the output: Fermi level = %.6f a.u.' % outinfo['fermi']
    if 'vbm' in outinfo:
        msg += ', VBM = %.6f a.u.' % outinfo['vbm']
    if 'cbm' in outinfo:
        msg += ', CBM = %.6f a.u.' % outinfo['cbm']
    if 'gap' in outinfo:
        msg += ', gap = %.3f eV' % (outinfo['gap'] * HARTREE2EV)
    ef_col2.info(msg + '. E = 0 at %.6f a.u.' % efermi)

elif ef_mode.startswith('Automatic'):
    if not cands:
        ef_col2.warning('No band gap could be detected (metallic spectrum?). '
                        'Upload the riper output above or use a custom value.')
        efermi = 0.0
    else:
        labels = ['band %d:  gap %.3f eV  (VBM at %.6f a.u.)'
                  % (i + 1, g, float(bmax[i])) for i, g, s in cands]
        pick = ef_col2.selectbox(
            'Detected gap candidates (flat core/semicore bands are already excluded; '
            'the first entry is the best guess)', labels, index=0)
        ivbm = cands[labels.index(pick)][0]
        nocc = ef_col2.number_input(
            'or set the number of occupied bands directly (0 = use the candidate above)',
            min_value=0, max_value=B['nbands'], value=0,
            help='For RKS this is the number of closed shells (= electrons per cell / 2), '
                 'see the $closed shells data group in control.')
        if nocc > 0:
            ivbm = min(max(int(nocc) - 1, 0), B['nbands'] - 2)
        efermi = float(bmax[ivbm])
        gap_ev = (bmin[ivbm + 1] - bmax[ivbm]) * HARTREE2EV
        ef_col2.info('E = 0 at the VBM: band %d, E = %.6f a.u.  →  gap %.3f eV.'
                     % (ivbm + 1, efermi, gap_ev))
        if gap_ev > 15.0 or gap_ev < 0.05:
            ef_col2.warning('The selected gap looks suspicious (%.2f eV). For metals or '
                            'unusual spectra the automatic detection cannot work - upload '
                            'the riper output or use a custom value.' % gap_ev)

else:
    efermi = ef_col2.number_input('Energy shift (a.u., e.g. the Fermi level printed by riper)',
                                  value=float(bmax[cands[0][0]]) if cands else 0.0,
                                  format='%.10f')
    ef_col2.caption('= %.4f eV. All bands are plotted relative to this value.'
                    % (efermi * HARTREE2EV))

# ----------------------------------------------------------------------------------------
# 2. What to plot
# ----------------------------------------------------------------------------------------
st.subheader('2. Bands and projections')

# default window: zoom on the gap region, not on the core bands
opt1, opt2, opt3, opt4 = st.columns(4)
emin = opt1.number_input('Energy window: min (eV)', value=-10.0, step=1.0)
emax = opt2.number_input('Energy window: max (eV)', value=10.0, step=1.0)
if B['nspin'] == 2:
    spin_pick = opt3.selectbox('Spin', ['alpha', 'beta', 'both'])
else:
    spin_pick = 'alpha'
    opt3.selectbox('Spin', ['closed shell'], disabled=True)
show_all = opt4.checkbox('Show all bands (no window filter)', value=False,
                         help='Draws every band including the core states. Can be slow and '
                              'usually not what you want.')
if emax <= emin:
    st.error('The energy window maximum must be larger than the minimum.')
    st.stop()

picked_buckets = []
if B['projected']:
    default_sel = [b for b in B['buckets'] if b.endswith('total') and b.startswith('element')]
    if not default_sel:
        default_sel = [b for b in B['buckets'] if b.endswith('total')][:4]
    picked_buckets = st.multiselect(
        'Projection buckets to overlay (marker area = Mulliken weight). Pick e.g. the '
        '*element ... s / p* buckets to see the orbital character of each element.',
        B['buckets'], default=default_sel[:4])

    with st.expander('How to interpret the fat-band markers', expanded=False):
        st.markdown("""
The **area of a circle is the Mulliken weight** of that bucket for that particular state
(band *n* at k-point *k*): half the area means half the weight. The size therefore changes
*along* a band whenever the character of the state changes with k — this k-dependence is
exactly the information the projection adds over a plain band plot.

**Example (MgO):** the lowest conduction band is a ~92 % Mg 3s state at Γ (large Mg
marker), but away from Γ it hybridizes — the Mg weight drops to ~82-84 % while the O weight
grows to ~16-18 %, and within Mg the character shifts from pure s to strongly p-mixed at M.
A shrinking marker of one element is always compensated by a growing marker of another.

Reading guide:

- **Sum rule.** The weights of one state over *all* atoms add up to 1. If your selected
  buckets cover the whole cell (e.g. the *total* buckets of all elements), their marker
  areas at any point represent fractions of one; a small marker of one bucket implies a
  large one of another. With a partial selection (a single atom, one l channel) the areas
  simply show that fraction of the state.
- **Overlapping circles.** Buckets are drawn in the order of the legend, semi-transparent,
  on top of the band line. A large circle with a small differently-colored dot at its
  center is a state dominated by the first bucket with a small admixture of the second -
  the layering is not a separate encoding.
- **Band crossings and degeneracies.** Where bands touch (degenerate multiplets, e.g. at
  high-symmetry points), the weights of the *individual* partner bands are arbitrary - only
  their sum over the touching bands is well defined. Erratic marker sizes right at a
  crossing are therefore not meaningful; a fraction of an eV away from the crossing they
  are.
- **Mulliken caveat.** Weights are Mulliken populations: basis-set dependent, and single
  shells can be slightly negative (such values are clipped to zero for drawing).
""")

# ----------------------------------------------------------------------------------------
# 3. Plot appearance (same spirit as the RT-TDDFT spectrum plotter)
# ----------------------------------------------------------------------------------------
st.subheader('3. Plot options')

t1, t2, t3 = st.columns(3)
plot_title = t1.text_input('Title', 'Band Structure')
ytitle = t2.text_input('y-axis label', 'E - E$_{VBM}$ (eV)' if ef_mode.startswith('Automatic')
                       else 'E - E$_F$ (eV)')
klabels_str = t3.text_input('High-symmetry point labels (comma separated, in path order)',
                            'Γ,X,M,Γ')

f1, f2, f3, f4 = st.columns(4)
figsize_w = f1.number_input('Figure width', value=9.0, step=0.5)
figsize_h = f2.number_input('Figure height', value=6.5, step=0.5)
band_color = f3.color_picker('Band line color', '#4477aa')
line_width = f4.slider('Band line width', 0.2, 5.0, 1.2, step=0.1)

g1, g2, g3, g4 = st.columns(4)
line_style = g1.selectbox('Band line style', ['solid', 'dashed', 'dotted', 'dashdot'])
marker_scale = g2.slider('Fat-band marker scale', 5, 200, 60,
                         disabled=not B['projected'])
marker_alpha = g3.slider('Fat-band marker opacity', 0.1, 1.0, 0.85, step=0.05,
                         disabled=not B['projected'])
transparency = g4.checkbox('Transparent background', value=False)

# ----------------------------------------------------------------------------------------
# Build the plot
# ----------------------------------------------------------------------------------------
kpath = B['kpath']
ticks = detect_ticks(B['kxyz'], kpath)
klabels = [s.strip() for s in klabels_str.split(',')] if klabels_str.strip() else []
if klabels and len(klabels) != len(ticks):
    st.warning('The path has %d segment boundaries but %d labels were given; the labels are '
               'skipped. (Coinciding segment ends count once.)' % (len(ticks), len(klabels)))
    klabels = []

spins = {'alpha': [1], 'beta': [2], 'both': [1, 2], 'closed shell': [1]}[spin_pick]
styles = {1: line_style, 2: 'dashed' if line_style == 'solid' else 'dotted'}

# Refuse to draw an empty frame: with a wrong energy reference the window may miss all bands
E0 = (B['E'][spins[0]] - efermi) * HARTREE2EV
if not show_all and not np.any((E0.min(axis=1) <= emax) & (E0.max(axis=1) >= emin)):
    enear = E0[np.argmin(np.abs(E0).min(axis=1))]
    st.error('No band lies inside the energy window %.1f ... %.1f eV - the plot would be '
             'empty. The nearest band spans %.2f ... %.2f eV. Check the energy reference '
             '(section 1) or widen the window.' % (emin, emax, enear.min(), enear.max()))
    st.stop()

fig, ax = plt.subplots(figsize=(figsize_w, figsize_h))

nshown = 0
for spin in spins:
    Eev = (B['E'][spin] - efermi) * HARTREE2EV
    if show_all:
        keep = np.arange(B['nbands'])
    else:
        keep = np.flatnonzero((Eev.min(axis=1) <= emax) & (Eev.max(axis=1) >= emin))
    nshown = max(nshown, len(keep))

    # all band lines in one LineCollection: fast even for very many bands
    segs = [np.column_stack([kpath, Eev[ib]]) for ib in keep]
    ax.add_collection(LineCollection(segs, colors=band_color, linewidths=line_width,
                                     linestyles=styles[spin], zorder=2))

    if B['projected'] and picked_buckets and spin == spins[0]:
        for ic, lab in enumerate(picked_buckets):
            iw = B['buckets'].index(lab)
            Wb = B['W'][spin][keep][:, :, iw]
            ww = np.clip(Wb, 0.0, None).ravel()
            mask = ww * marker_scale > 0.5          # skip invisible markers
            xx = np.tile(kpath, len(keep))[mask]
            yy = Eev[keep].ravel()[mask]
            ax.scatter(xx, yy, s=ww[mask] * marker_scale, color=PALETTE[ic % len(PALETTE)],
                       alpha=marker_alpha, edgecolors='none', zorder=3 + ic, label=lab)

nskip = B['nbands'] - nshown
if nskip > 0:
    st.caption('%d of %d bands are inside the energy window; %d bands (mostly core states) '
               'are not drawn. Widen the window or tick *Show all bands* to see them.'
               % (nshown, B['nbands'], nskip))

for t in ticks:
    ax.axvline(t, color='gray', lw=0.8, zorder=1)
ax.axhline(0.0, color='black', ls='--', lw=1.2, zorder=1)
if klabels:
    ax.set_xticks(ticks)
    ax.set_xticklabels(klabels, fontsize=15)
else:
    ax.set_xticks(ticks)
    ax.set_xticklabels([''] * len(ticks))
ax.tick_params(axis='y', labelsize=13)
ax.set_xlim(kpath[0], kpath[-1])
ax.set_ylim(emin, emax)
ax.set_ylabel(ytitle, fontsize=15)
ax.set_title(plot_title, fontsize=17)
if B['projected'] and picked_buckets:
    ax.legend(fontsize=10, markerscale=1.0, loc='best')
if spin_pick == 'both':
    ax.plot([], [], color=band_color, ls=styles[1], label='alpha')
    ax.plot([], [], color=band_color, ls=styles[2], label='beta')
    ax.legend(fontsize=10, loc='best')

fig.tight_layout()
st.pyplot(fig)

# ----------------------------------------------------------------------------------------
# Downloads
# ----------------------------------------------------------------------------------------
buf = io.BytesIO()
fig.savefig(buf, format='png', dpi=300, transparent=transparency)
d1, d2 = st.columns(2)
d1.download_button('Download plot as PNG file', buf.getvalue(), file_name='bands_plot.png',
                   mime='image/png')

# plotted data as CSV (bands inside the window, first selected spin)
Eev = (B['E'][spins[0]] - efermi) * HARTREE2EV
keep = np.arange(B['nbands']) if show_all else \
    np.flatnonzero((Eev.min(axis=1) <= emax) & (Eev.max(axis=1) >= emin))
out = pd.DataFrame({'|k|': kpath})
for ib in keep:
    out['band_%d_eV' % (ib + 1)] = Eev[ib]
d2.download_button('Download the plotted bands as CSV', out.to_csv(index=False).encode(),
                   file_name='bands_window.csv', mime='text/csv')
