"""RIPER band structure plotter, adapted from RIPER-Tools-for-TURBOMOLE
band_plotting_examples/3D_Si/bands_plotter_RKS.py (same parsing and plot style).

Two modes:
  normal    - plot bands.xyz exactly like the original script
  projected - plot bands_proj.xyz: same bands, plus circles whose AREA is the
              Mulliken weight of two chosen buckets (fat bands)

Usage:
  python bands_plotter_adapted.py normal    <bands.xyz>      <shift_au> <nk> <emin> <emax> <title> <out.png> <labels>
  python bands_plotter_adapted.py projected <bands_proj.xyz> <shift_au> <nk> <emin> <emax> <title> <out.png> <labels> <bucketA> <bucketB>
"""
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

######### INPUT STUFF (from the command line instead of edited constants)
mode = sys.argv[1]
filename = sys.argv[2]
energy_shift = float(sys.argv[3])          # User-provided energy shift value (in au)
num_k_points = int(sys.argv[4])            # Num k-points for each band
energy_range_min = float(sys.argv[5])      # Energy range for plotting (in eV)
energy_range_max = float(sys.argv[6])
plot_title = sys.argv[7]
outfile = sys.argv[8]
kpoint_labels_plot = sys.argv[9].split(',')
wanted = sys.argv[10:12] if mode == 'projected' else []

# Initialize variables to store data
up_spin_data = []
bucket_labels = {}
is_up_spin = True  # Flag to track whether we're reading up spin or down spin data

# Read the file
with open(filename, 'r') as file:
    for line in file:
        line = line.strip()

        if line.startswith('# bucket'):
            num, lab = line[8:].split(':', 1)
            bucket_labels[int(num)] = lab.strip()
            continue
        if line.startswith('#') or line == '':
            continue
        if line.startswith('closed/alpha shells:'):
            continue
        if line.startswith('beta shells:'):
            is_up_spin = False
            continue

        # Split the line into values
        values = line.split()

        if is_up_spin and len(values) >= 5:  # kx, ky, kz, |k|, energy (+ weights if projected)
            up_spin_data.append(values)

# Create dataframes
ncols = len(up_spin_data[0])
columns = ['kx', 'ky', 'kz', '|k|', 'energy'] + ['w%d' % i for i in range(1, ncols - 4)]
up_spin_df = pd.DataFrame(up_spin_data, columns=columns)
up_spin_df = up_spin_df.astype({c: float for c in columns})

# Calculate the number of bands for each spin
num_bands_up_spin = up_spin_df.shape[0] // num_k_points
print("Number of Bands :", num_bands_up_spin)

# Shift energies by the user-provided value and convert to eV
conversion_factor = 27.211324570273
up_spin_df['energy'] -= energy_shift
up_spin_df['energy_ev'] = up_spin_df['energy'] * conversion_factor

# k-point ticks: the path column of one band; a segment boundary is where the direction of the
# k vector changes (segment lengths may coincide, so spacing alone is not enough)
kpath = up_spin_df['|k|'][0:num_k_points].tolist()
kx = up_spin_df['kx'][0:num_k_points].tolist()
ky = up_spin_df['ky'][0:num_k_points].tolist()
kz = up_spin_df['kz'][0:num_k_points].tolist()
kpoint_ticks = [kpath[0]]
for j in range(1, num_k_points - 1):
    d1 = (kx[j] - kx[j - 1], ky[j] - ky[j - 1], kz[j] - kz[j - 1])
    d2 = (kx[j + 1] - kx[j], ky[j + 1] - ky[j], kz[j + 1] - kz[j])
    n1 = max((d1[0]**2 + d1[1]**2 + d1[2]**2)**0.5, 1e-30)
    n2 = max((d2[0]**2 + d2[1]**2 + d2[2]**2)**0.5, 1e-30)
    cosang = (d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2]) / (n1 * n2)
    if cosang < 1.0 - 1e-8:
        kpoint_ticks.append(kpath[j])
kpoint_ticks.append(kpath[-1])

# Find bands within the specified energy range
bands_to_plot = []
for band_index in range(num_bands_up_spin):
    band_data = up_spin_df[band_index * num_k_points: (band_index + 1) * num_k_points]
    if any((energy_range_min <= energy <= energy_range_max) for energy in band_data['energy_ev']):
        bands_to_plot.append(band_data)

plt.figure(figsize=(8.0, 6.0))

# Plot bands
for band_data in bands_to_plot:
    plt.plot(kpath, band_data['energy_ev'].tolist(), color='steelblue', linestyle='--', lw=2.0)

# Fat-band overlay: circle area = Mulliken weight of the chosen buckets
if mode == 'projected':
    fat_colors = ['#fc8d62', '#66c2a5']
    picks = []
    for w in wanted:
        hits = [i for i, lab in bucket_labels.items() if w in lab]
        if not hits:
            raise SystemExit('no bucket matches %r; available: %s' % (w, list(bucket_labels.values())))
        picks.append(hits[0])
    for ic, ib in enumerate(picks):
        col = 'w%d' % ib
        shown = False
        for band_data in bands_to_plot:
            ww = band_data[col].clip(lower=0.0)
            plt.scatter(kpath, band_data['energy_ev'], s=60.0 * ww, color=fat_colors[ic],
                        alpha=0.9, edgecolors='none', zorder=3 + ic,
                        label=None if shown else bucket_labels[ib])
            shown = True
    plt.legend(fontsize=14, loc='upper right', markerscale=1.3)

# Set custom k-point labels and ticks
plt.xticks(kpoint_ticks, kpoint_labels_plot, fontsize=16)
plt.yticks(fontsize=16)

plt.ylabel('Energy (eV)', fontsize=18)
plt.title(plot_title, fontsize=20)
plt.axhline(y=0, color='black', linestyle='--', linewidth=1.7)  # Add a horizontal line at y=0
# Disable horizontal grid lines
plt.grid(which='major', axis='x', linestyle='-', linewidth=1.8, color='gray')
plt.xlim(kpoint_ticks[0], kpoint_ticks[-1])
plt.ylim(energy_range_min, energy_range_max)
plt.tight_layout()
plt.savefig(outfile, dpi=160)
print('wrote', outfile)
