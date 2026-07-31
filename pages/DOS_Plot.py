import io
import re
import zipfile

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Set page config
st.set_page_config(page_title='Density of States Plotter for RIPER', layout='wide', page_icon="⚛️",
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
PALETTE = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494',
           '#b3b3b3', '#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02']
SPIN_SUFFIXES = ['_a+b', '_a-b', '_alpha', '_beta']
SPIN_LABEL = {'a+b': 'alpha + beta', 'a-b': 'alpha - beta', 'alpha': 'alpha',
              'beta': 'beta', 'none': 'spin restricted'}

# Files a user is likely to upload by mistake, and what to tell them
LOOKS_LIKE = [
    ('$coord', 'this is a **coord** file (the structure), not a DOS file'),
    ('$cell', 'this is a **coord** file (the structure), not a DOS file'),
    ('$dosper', 'this is the **control** file. It *requests* the DOS; the DOS itself is '
                'written to separate files when you run `riper`'),
    ('$pdosper', 'this is the **control** file. It *requests* the PDOS; the PDOS itself is '
                 'written to separate files when you run `riper`'),
    ('$scfmo', 'this is an **mos**/orbital file, not a DOS file'),
    ('closed/alpha shells:', 'this is **bands.xyz** (a band structure), not a DOS file. '
                             'Use the *Band Structure* page for it'),
    ('$grad', 'this is a **gradient** file, not a DOS file'),
]


# ----------------------------------------------------------------------------------------
# Parsing
# ----------------------------------------------------------------------------------------
def _classify_header(line):
    """Return (kind, label, sortkey) from a RIPER DOS/PDOS identification line."""
    m = re.search(r'PDOS\s+shell\s+(\d+)\s*\((\w)\)\s+of\s+atom\s+(\d+)\s+(\S+)', line)
    if m:
        sh, lq, iat, sym = int(m.group(1)), m.group(2), int(m.group(3)), m.group(4)
        return 'shell', '%s atom %d, shell %d (%s)' % (sym.capitalize(), iat, sh, lq), (3, iat, sh)
    m = re.search(r'PDOS\s+element\s+(\S+)', line)
    if m:
        sym = m.group(1)
        return 'element', '%s (element)' % sym.capitalize(), (1, 0, sym)
    m = re.search(r'PDOS\s+atom\s+(\d+)\s+(\S+)', line)
    if m:
        iat, sym = int(m.group(1)), m.group(2)
        return 'atom', '%s atom %d' % (sym.capitalize(), iat), (2, iat, 0)
    if 'DOS' in line:
        return 'total', 'Total DOS', (0, 0, 0)
    return None, None, (9, 0, 0)


def _spin_from_name(name):
    base = name.strip()
    for suf in SPIN_SUFFIXES:
        if base.endswith(suf):
            return suf[1:], base[:-len(suf)]
    return 'none', base


def parse_riper_dos(text, filename=''):
    """Parse a RIPER dos/pdos file.

    Returns a dict with the parsed content, or raises ValueError with a message that
    tells the user what to do about it.
    """
    if not text or not text.strip():
        raise ValueError('the file is empty')

    lines = text.splitlines()

    # Catch the common "wrong file" mistakes before anything else
    head = '\n'.join(lines[:60])
    for token, hint in LOOKS_LIKE:
        if token in head:
            raise ValueError(hint)

    kind, label, sortkey = None, None, (9, 0, 0)
    columns, ef, shifted = None, None, False
    rows, bad_rows = [], []

    for iline, raw in enumerate(lines, start=1):
        s = raw.strip()
        if not s:
            continue
        if s.startswith('#'):
            body = s.lstrip('#').strip()
            if kind is None:
                k, lab, sk = _classify_header(body)
                if k is not None:
                    kind, label, sortkey = k, lab, sk
                    continue
            if 'total' in body.split():
                tok = body.split()
                if len(tok) >= 2 and tok[0] == 'E' and tok[1] == '-':
                    shifted, columns = True, tok[3:]
                else:
                    columns = tok[1:]
                continue
            m = re.search(r'Ef[^=]*=\s*(-?\d+\.?\d*(?:[eEdD][+-]?\d+)?)', body)
            if m and ef is None:
                ef = float(m.group(1).replace('D', 'E').replace('d', 'e'))
            continue

        parts = s.split()
        try:
            vals = [float(p.replace('D', 'E').replace('d', 'e')) for p in parts]
        except ValueError:
            bad_rows.append(iline)
            continue
        rows.append(vals)

    if not rows:
        if kind is None and columns is None:
            raise ValueError('this does not look like a RIPER DOS file. Expected a `#` header '
                             'such as `#     DOS` or `#     PDOS  element o  (Mulliken)` '
                             'followed by columns of numbers')
        raise ValueError('the header was found but the file contains no data rows. It may be '
                         'truncated, or the calculation may not have finished writing it')

    width = max(len(r) for r in rows)
    rows = [r for r in rows if len(r) == width]
    data = np.array(rows, dtype=float)

    if columns is None:
        # Header missing (e.g. the user stripped the comment lines when copying)
        columns = ['total'] + ['col%d' % i for i in range(2, width)]
        headerless = True
    else:
        headerless = False
        if len(columns) != width - 1:
            # Trust the data, pad/trim the names
            if len(columns) < width - 1:
                columns = columns + ['col%d' % i for i in range(len(columns) + 2, width + 1)]
            else:
                columns = columns[:width - 1]

    if kind is None:
        kind, label, sortkey = 'total', 'Total DOS', (0, 0, 0)
        if filename:
            kind2, label2, sortkey2 = _classify_header(filename.replace('_', ' '))
            if kind2 and kind2 != 'total':
                kind, label, sortkey = kind2, label2, sortkey2

    spin, base = _spin_from_name(filename) if filename else ('none', '')

    return {
        'kind': kind, 'label': label, 'sortkey': sortkey, 'filename': filename or '(pasted)',
        'basename': base, 'spin': spin, 'columns': columns, 'shifted': shifted, 'ef': ef,
        'energy': data[:, 0], 'values': data[:, 1:], 'npts': data.shape[0],
        'bad_rows': bad_rows, 'headerless': headerless,
    }


def grid_signature(rec):
    return (rec['npts'], round(float(rec['energy'][0]), 8), round(float(rec['energy'][-1]), 8))


# ----------------------------------------------------------------------------------------
# Page
# ----------------------------------------------------------------------------------------
st.title("Density of States Plotter")
st.write('Plots the **total DOS** and the **projected DOS (PDOS)** written by `riper`. '
         'Upload as many files as you like — atoms, elements and shells are recognised '
         'automatically from the file headers.')

with st.expander("How to produce these files with RIPER", expanded=False):
    st.write('**Total DOS only** &ndash; add to your `control` file:')
    st.code("$dosper width=0.01 npt=3000 fermishift", language='text')
    st.write('**Projected DOS** &ndash; add `$pdosper` as well. It reuses the energy grid of '
             '`$dosper`, so both keywords are needed:')
    st.code("""$dosper  width=0.01 npt=3000 fermishift
$pdosper
  atoms    1,3-5        # optional, default: all atoms
  elements o,ti         # optional, default: all elements
  resolve  l            # atom | l | shell   (default: l)""", language='text')
    st.write('Then run `riper` (a normal single point run, or `-proper` to reuse converged orbitals):')
    st.code('nohup riper -proper > dos.out &', language='shell')
    st.write('This writes `dos` (or `dos_a+b`, `dos_a-b`, `dos_alpha`, `dos_beta` for UKS) and, '
             'with `$pdosper`, files named `pdos_at<i>_<sym>`, `pdos_el_<sym>` and '
             '`pdos_at<i>_<sym>_sh<m>`. **Upload all of them at once.**')

mode = st.radio("How do you want to provide the data?",
                ["Upload files (recommended, works for total DOS and PDOS)",
                 "Paste the contents of a single file"],
                index=0)

records, problems = [], []

if mode.startswith("Upload"):
    ups = st.file_uploader(
        "Drop your `dos` / `pdos_*` files here (you can select many at once, or a .zip of them)",
        accept_multiple_files=True)

    raw_items = []
    for up in ups or []:
        try:
            blob = up.read()
        except Exception as exc:                                    # pragma: no cover
            problems.append((up.name, 'could not be read (%s)' % exc))
            continue
        if up.name.lower().endswith('.zip'):
            try:
                with zipfile.ZipFile(io.BytesIO(blob)) as zf:
                    for member in zf.namelist():
                        if member.endswith('/'):
                            continue
                        raw_items.append((member.split('/')[-1], zf.read(member)))
            except Exception as exc:
                problems.append((up.name, 'is not a readable zip archive (%s)' % exc))
            continue
        raw_items.append((up.name, blob))

    seen = {}
    for name, blob in raw_items:
        try:
            text = blob.decode('utf-8', errors='replace')
        except Exception:
            problems.append((name, 'is not a text file'))
            continue
        if name in seen:
            continue
        seen[name] = True
        try:
            records.append(parse_riper_dos(text, name))
        except ValueError as exc:
            problems.append((name, str(exc)))

else:
    st.write('Paste one file below. For several files use the upload mode instead &ndash; '
             'the file *names* are what identify the spin channel.')
    pasted_name = st.text_input("File name (optional, helps identify the file)", value="")
    pasted = st.text_area("File contents", height=300)
    if pasted.strip():
        try:
            records.append(parse_riper_dos(pasted, pasted_name.strip()))
        except ValueError as exc:
            problems.append((pasted_name.strip() or 'pasted text', str(exc)))


# ----------------------------------------------------------------------------------------
# Report what could not be used
# ----------------------------------------------------------------------------------------
for name, msg in problems:
    st.error("**%s** could not be used: %s." % (name, msg))

if problems and not records:
    st.info("Nothing could be plotted yet. The files written by `riper` start with a comment "
            "line such as `#     DOS` or `#     PDOS  element o  (Mulliken)`, followed by a "
            "column header and then two or more columns of numbers.")

if not records:
    if not problems:
        if mode.startswith("Upload"):
            st.info("⬆️ Upload at least one `dos` or `pdos_*` file to get started. Open "
                    "*How to produce these files with RIPER* above if you do not have them yet.")
        else:
            st.info("⬆️ Paste the contents of one `dos` or `pdos_*` file above, then press "
                    "**Ctrl+Enter** (⌘+Enter on a Mac) to apply it. Open *How to produce these "
                    "files with RIPER* above if you do not have such a file yet.")
    st.stop()


# ----------------------------------------------------------------------------------------
# Consistency checks
# ----------------------------------------------------------------------------------------
for rec in records:
    if rec['bad_rows']:
        st.warning("**%s**: skipped %d line(s) that were not numeric (first at line %d)."
                   % (rec['filename'], len(rec['bad_rows']), rec['bad_rows'][0]))
    if rec['headerless']:
        st.warning("**%s**: no `#` header found, so the columns could not be named. The first "
                   "column is assumed to be the energy. Upload the original file to get proper "
                   "labels." % rec['filename'])

grids = {}
for rec in records:
    grids.setdefault(grid_signature(rec), []).append(rec['filename'])

if len(grids) > 1:
    st.error("The uploaded files do not share the same energy grid, so they cannot be compared "
             "in one plot. All files must come from **one** `riper` run with **one** `$dosper` "
             "setting.")
    for (npts, e0, e1), names in grids.items():
        st.write("- %d points from %.5f to %.5f: `%s`" % (npts, e0, e1, '`, `'.join(names)))
    st.info("Re-run `riper -proper` once with the `$dosper` settings you want, then upload the "
            "files it produces together.")
    st.stop()

shift_states = set(r['shifted'] for r in records)
if len(shift_states) > 1:
    st.warning("Some files are referenced to the Fermi level and some are not. They were written "
               "by different runs (`fermishift` vs `nofermishift`); the energy axes do not match.")

spins_present = sorted(set(r['spin'] for r in records))
if 'none' in spins_present and len(spins_present) > 1:
    st.warning("You mixed spin-restricted files (no suffix) with spin-resolved ones "
               "(`_alpha`, `_beta`, `_a+b`, `_a-b`). Pick one channel below.")


# ----------------------------------------------------------------------------------------
# Spin channel selection
# ----------------------------------------------------------------------------------------
st.subheader("1. Data that was recognised")

summary = pd.DataFrame([{
    'file': r['filename'],
    'type': {'total': 'total DOS', 'element': 'element', 'atom': 'atom', 'shell': 'shell'}[r['kind']],
    'label': r['label'],
    'spin': SPIN_LABEL.get(r['spin'], r['spin']),
    'columns': ', '.join(r['columns']),
    'points': r['npts'],
} for r in sorted(records, key=lambda r: (r['spin'], r['sortkey']))])
st.dataframe(summary, use_container_width=True, hide_index=True)

if len(spins_present) > 1:
    nice = [SPIN_LABEL.get(s, s) for s in spins_present]
    pick = st.selectbox("Spin channel to plot", nice,
                        index=nice.index('alpha') if 'alpha' in nice else 0)
    if pick == 'alpha' and 'beta' in spins_present:
        mirror = st.checkbox("Plot beta mirrored below the axis (spin-up / spin-down plot)",
                             value=True)
    else:
        mirror = False
    chosen = spins_present[nice.index(pick)]
else:
    chosen, mirror = spins_present[0], False

active = [r for r in records if r['spin'] == chosen]
beta = [r for r in records if r['spin'] == 'beta'] if mirror else []

ef_vals = [r['ef'] for r in active if r['ef'] is not None]
ef = ef_vals[0] if ef_vals else None
is_shifted = all(r['shifted'] for r in active)


# ----------------------------------------------------------------------------------------
# Curve selection
# ----------------------------------------------------------------------------------------
st.subheader("2. Curves to plot")

curves = []
for rec in sorted(active, key=lambda r: r['sortkey']):
    for ic, cname in enumerate(rec['columns']):
        curves.append({'rec': rec, 'ic': ic,
                       'name': '%s — %s' % (rec['label'], cname),
                       'is_total_col': (cname == 'total'),
                       'kind': rec['kind']})

names = [c['name'] for c in curves]

presets = {
    "Elements (total)": [c['name'] for c in curves if c['is_total_col'] and c['kind'] == 'element'],
    "Atoms (total)": [c['name'] for c in curves if c['is_total_col'] and c['kind'] == 'atom'],
    "Shells": [c['name'] for c in curves if c['kind'] == 'shell'],
    "Projections split by l": [c['name'] for c in curves
                               if not c['is_total_col'] and c['kind'] != 'total'],
    "Total DOS split by l": [c['name'] for c in curves
                             if not c['is_total_col'] and c['kind'] == 'total'],
}
presets = {k: v for k, v in presets.items() if v}
if not presets:
    # Only bare total columns are available (e.g. a single pasted shell file)
    presets = {"All curves": [c['name'] for c in curves]}
presets["Custom"] = []

pnames = list(presets.keys())
preset = st.radio("Quick selection", pnames, index=0, horizontal=True,
                  help="Pick a starting set, then add or remove individual curves below.")

# Re-key the multiselect on the preset so switching preset actually resets the widget
picked = st.multiselect("Curves", names, default=presets[preset][:14],
                        key="curves_%s_%s" % (preset, chosen))
st.caption("%d curve(s) available in total." % len(names))

if not picked:
    st.info("Select at least one curve above to draw the plot. The **— total** entries are the "
            "summed contribution of an atom or element; the **— s / p / d** entries split it by "
            "angular momentum.")
    st.stop()

total_rec = next((r for r in active if r['kind'] == 'total'), None)


# ----------------------------------------------------------------------------------------
# Sum rule check
# ----------------------------------------------------------------------------------------
if total_rec is not None:
    els = [r for r in active if r['kind'] == 'element']
    ats = [r for r in active if r['kind'] == 'atom']
    check_set, check_what = (els, 'elements') if els else ((ats, 'atoms') if ats else (None, None))
    if check_set:
        ssum = np.zeros_like(total_rec['values'][:, 0])
        for r in check_set:
            ssum = ssum + r['values'][:, 0]
        resid = float(np.max(np.abs(ssum - total_rec['values'][:, 0])))
        peak = float(np.max(np.abs(total_rec['values'][:, 0]))) or 1.0
        if resid <= 5e-5 * max(1.0, peak):
            st.success("Sum rule satisfied: the %s add up to the total DOS "
                       "(largest deviation %.2e, i.e. the printing precision of the files)."
                       % (check_what, resid))
        else:
            st.warning("The %s do not add up to the total DOS (largest deviation %.3g). "
                       "This normally means some files are missing from the upload &ndash; the "
                       "sum only closes when **every** atom or **every** element of the cell is "
                       "included." % (check_what, resid))


# ----------------------------------------------------------------------------------------
# Plot options
# ----------------------------------------------------------------------------------------
st.subheader("3. Plot")

o1, o2, o3, o4 = st.columns(4)
unit = o1.selectbox("Energy unit", ["eV", "Hartree"], index=0)
show_total = o2.checkbox("Shade the total DOS", value=total_rec is not None,
                         disabled=total_rec is None,
                         help=None if total_rec is not None else
                         "Upload the `dos` file as well to enable this")
line_width = o3.slider("Line width", 0.4, 6.0, 2.1, step=0.1)
transparency = o4.checkbox("Transparent background", value=False)

p1, p2, p3, p4 = st.columns(4)
fig_x = p1.number_input("Figure width", value=10.0, step=0.5)
fig_y = p2.number_input("Figure height", value=5.0, step=0.5)
is_grid = p3.checkbox("Grid", value=False)
draw_ef = p4.checkbox("Line at the Fermi level", value=True)

conv = HARTREE2EV if unit == "eV" else 1.0


def energy_of(rec):
    e = rec['energy'] * conv
    if not rec['shifted'] and rec['ef'] is not None:
        e = e - rec['ef'] * conv
    return e


def dos_of(rec, ic):
    return rec['values'][:, ic] / conv


emin_d = float(min(energy_of(r).min() for r in active))
emax_d = float(max(energy_of(r).max() for r in active))
r1, r2 = st.columns(2)
xmin = r1.number_input("Energy axis: minimum", value=round(emin_d, 3))
xmax = r2.number_input("Energy axis: maximum", value=round(emax_d, 3))
if xmax <= xmin:
    st.error("The maximum of the energy axis must be larger than the minimum. "
             "Reset them to %.3f and %.3f to see the whole range." % (emin_d, emax_d))
    st.stop()

t1, t2 = st.columns(2)
title = t1.text_input("Title", "Density of States")
ylabel = t2.text_input("y-axis label", "DOS (states/%s/cell)" % unit)

fig, ax = plt.subplots(figsize=(fig_x, fig_y))
ax.grid(is_grid)

if show_total and total_rec is not None:
    ax.fill_between(energy_of(total_rec), dos_of(total_rec, 0), color='lightgrey',
                    label='total DOS', zorder=0)
    if mirror:
        btot = next((r for r in beta if r['kind'] == 'total'), None)
        if btot is not None:
            ax.fill_between(energy_of(btot), -dos_of(btot, 0), color='lightgrey', zorder=0)

for i, cname in enumerate(picked):
    c = next(cc for cc in curves if cc['name'] == cname)
    ax.plot(energy_of(c['rec']), dos_of(c['rec'], c['ic']),
            color=PALETTE[i % len(PALETTE)], lw=line_width, label=cname, zorder=2)
    if mirror:
        bmatch = next((r for r in beta if r['basename'] == c['rec']['basename']), None)
        if bmatch is not None and c['ic'] < bmatch['values'].shape[1]:
            ax.plot(energy_of(bmatch), -dos_of(bmatch, c['ic']),
                    color=PALETTE[i % len(PALETTE)], lw=line_width, zorder=2)

if draw_ef:
    ax.axvline(0.0 if is_shifted else (ef * conv if ef is not None else 0.0),
               color='black', linestyle='--', lw=0.6)
if mirror:
    ax.axhline(0.0, color='black', lw=0.6)

ax.set_xlabel(('E - E$_F$ (%s)' if is_shifted else 'Energy (%s)') % unit)
ax.set_ylabel(ylabel)
ax.set_xlim(xmin, xmax)
if not mirror:
    ax.set_ylim(bottom=0.0)
ax.set_title(title)
ncol = 2 if len(picked) > 8 else 1
ax.legend(fontsize=8, ncol=ncol)
fig.tight_layout()

st.pyplot(fig)

if not is_shifted and ef is None:
    st.info("These files were written with `nofermishift` and carry no Fermi level, so the "
            "energy axis is absolute. Add `fermishift` to `$dosper` to reference the "
            "energies to E$_F$.")


# ----------------------------------------------------------------------------------------
# Downloads
# ----------------------------------------------------------------------------------------
st.subheader("4. Download")

buf = io.BytesIO()
fig.savefig(buf, format='png', dpi=400, transparent=transparency)
d1, d2 = st.columns(2)
d1.download_button('Download plot as PNG', buf.getvalue(), file_name='dos_plot.png',
                   mime='image/png')

out = pd.DataFrame({'energy_%s' % unit: energy_of(active[0])})
for cname in picked:
    c = next(cc for cc in curves if cc['name'] == cname)
    out[cname] = dos_of(c['rec'], c['ic'])
d2.download_button('Download the plotted data as CSV', out.to_csv(index=False).encode(),
                   file_name='dos_data.csv', mime='text/csv')

with st.expander("Show the plotted data as a table"):
    st.dataframe(out, use_container_width=True)
