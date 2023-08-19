import pandas as pd
import matplotlib.pyplot as plt


######### INPUT STUFF
# Define the filename
filename = "C:/Users/manas/Documents/LK99/github/LK99/TURBOMOLE_DFT_Calculations/bands_TPSS.xyz"

# Define the energy range for plotting (in eV)
energy_range_min = -1  # Modify this according to your desired range
energy_range_max = 0.1   # Modify this according to your desired range
# User-provided energy shift value
energy_shift =  -0.088521#-0.094143#-0.035363#-0.09
# Num k-points for each band
num_k_points = 360
# K-point labels
kpoint_labels_plot = ['Γ', 'M', 'K', 'Γ', 'A', 'L', 'H', 'A|L', 'M|K', 'H']
kpoint_ticks = [0.0000000, 0.1944399225, 0.3709054431, 0.7182621428, 0.9516852395, 1.1518848865, 1.3283504070, 1.6757071068, 1.9091302034, 2.1395191840]



# Initialize variables to store data
up_spin_data = []
down_spin_data = []
is_up_spin = True  # Flag to track whether we're reading up spin or down spin data

# Read the file
with open(filename, 'r') as file:
    for line in file:
        line = line.strip()
        
        if line == 'beta shells:':
            is_up_spin = False
            continue
        
        if line.startswith('closed/alpha shells:') or line == '':
            continue
        
        # Split the line into values
        values = line.split()
        
        if len(values) == 5:  # Assuming each line has kx, ky, kz, |k|, and energy
            if is_up_spin:
                up_spin_data.append(values)
            else:
                down_spin_data.append(values)

# Create dataframes
up_spin_df = pd.DataFrame(up_spin_data, columns=['kx', 'ky', 'kz', '|k|', 'energy'])
down_spin_df = pd.DataFrame(down_spin_data, columns=['kx', 'ky', 'kz', '|k|', 'energy'])

# Convert columns to appropriate data types
up_spin_df = up_spin_df.astype({'kx': float, 'ky': float, 'kz': float, '|k|': float, 'energy': float})
down_spin_df = down_spin_df.astype({'kx': float, 'ky': float, 'kz': float, '|k|': float, 'energy': float})

# Print the dataframes
print("Up Spin Data:")
print(up_spin_df)

print("\nDown Spin Data:")
print(down_spin_df)

# Calculate the number of bands for each spin
num_bands_up_spin = up_spin_df.shape[0] // num_k_points
num_bands_down_spin = down_spin_df.shape[0] // num_k_points

print("Number of Bands - Up Spin:", num_bands_up_spin)
print("Number of Bands - Down Spin:", num_bands_down_spin)



# Shift energies by the user-provided value
up_spin_df['energy'] -= energy_shift
down_spin_df['energy'] -= energy_shift

# Convert energies to eV
conversion_factor = 27.211324570273
up_spin_df['energy_ev'] = up_spin_df['energy'] * conversion_factor
down_spin_df['energy_ev'] = down_spin_df['energy'] * conversion_factor

# Print the modified dataframes
print("\nUp Spin Data (Energy Shifted and in eV):")
print(up_spin_df)

print("\nDown Spin Data (Energy Shifted and in eV):")
print(down_spin_df)


# Find bands within the specified energy range for both spins
up_spin_bands_to_plot = []
down_spin_bands_to_plot = []

for band_index in range(num_bands_up_spin):
    band_data = up_spin_df[band_index * num_k_points : (band_index + 1) * num_k_points]
    if any((energy_range_min <= energy <= energy_range_max) for energy in band_data['energy_ev']):
        up_spin_bands_to_plot.append(band_data['energy_ev'].tolist())


for band_index in range(num_bands_down_spin):
    band_data = down_spin_df[band_index * num_k_points : (band_index + 1) * num_k_points]
    if any((energy_range_min <= energy <= energy_range_max) for energy in band_data['energy_ev']):
        down_spin_bands_to_plot.append(band_data['energy_ev'].tolist())


# Plot bands for up spin
for band_data in up_spin_bands_to_plot:
    plt.plot(up_spin_df['|k|'][0:num_k_points], band_data, label=f'Band {band_index + 1}', color='steelblue', linestyle='--', lw=2.0)

# Plot bands for up spin
for band_data in down_spin_bands_to_plot:
    plt.plot(down_spin_df['|k|'][0:num_k_points], band_data, label=f'Band {band_index + 1}', color='orange', lw=2.0)


# Set custom k-point labels and ticks
plt.xticks(kpoint_ticks, kpoint_labels_plot, fontsize=16)


plt.yticks(fontsize=16)

# plt.xlabel('|k|', fontsize=18)
plt.ylabel('Energy (eV)', fontsize=18)
plt.title('Band Structure (TPSS)', fontsize=20)
plt.axhline(y=0, color='black', linestyle='--', linewidth=1.7)  # Add a horizontal line at y=0
# plt.legend()
# Disable horizontal grid lines
plt.grid(which='major', axis='x', linestyle='-', linewidth=1.8, color='gray')
plt.xlim(kpoint_ticks[0], kpoint_ticks[-1])
plt.tight_layout()
plt.show()
