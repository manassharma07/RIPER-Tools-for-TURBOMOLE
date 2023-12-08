import pandas as pd
import matplotlib.pyplot as plt


######### INPUT STUFF
# Define the filename
filename = "bands.xyz"

# Define the energy range for plotting (in eV)
energy_range_min = -2.5  # Modify this according to your desired range
energy_range_max = 3   # Modify this according to your desired range
# User-provided energy shift value (in au)
energy_shift =  -0.237887
# Num k-points for each band
num_k_points = 400
# K-point labels
kpoint_labels_plot = ['Γ', 'X', 'W', 'K', 'Γ', 'L', 'U', 'W', 'L', 'K|U', 'X']
kpoint_ticks = [0.0000000, 0.5955129971, 0.9085390597, 1.1267194809, 1.7637550049, 2.2956802070, 2.6735797819, 2.8934760813, 3.3199659659, 3.6954389227, 3.9153352221]
plot_title = 'Band Structure (PBE)'



# Initialize variables to store data
up_spin_data = []
is_up_spin = True  # Flag to track whether we're reading up spin or down spin data

# Read the file
with open(filename, 'r') as file:
    for line in file:
        line = line.strip()
        
        
        if line.startswith('closed/alpha shells:') or line == '':
            continue
        
        # Split the line into values
        values = line.split()
        
        if len(values) == 5:  # Assuming each line has kx, ky, kz, |k|, and energy
            if is_up_spin:
                up_spin_data.append(values)
            
# Create dataframes
up_spin_df = pd.DataFrame(up_spin_data, columns=['kx', 'ky', 'kz', '|k|', 'energy'])

# Convert columns to appropriate data types
up_spin_df = up_spin_df.astype({'kx': float, 'ky': float, 'kz': float, '|k|': float, 'energy': float})

# Print the dataframes
print("Up Spin Data:")
print(up_spin_df)


# Calculate the number of bands for each spin
num_bands_up_spin = up_spin_df.shape[0] // num_k_points

print("Number of Bands :", num_bands_up_spin)



# Shift energies by the user-provided value
up_spin_df['energy'] -= energy_shift

# Convert energies to eV
conversion_factor = 27.211324570273
up_spin_df['energy_ev'] = up_spin_df['energy'] * conversion_factor

# Print the modified dataframes
print("\nUp Spin Data (Energy Shifted and in eV):")
print(up_spin_df)


# Find bands within the specified energy range 
up_spin_bands_to_plot = []
down_spin_bands_to_plot = []

for band_index in range(num_bands_up_spin):
    band_data = up_spin_df[band_index * num_k_points : (band_index + 1) * num_k_points]
    if any((energy_range_min <= energy <= energy_range_max) for energy in band_data['energy_ev']):
        up_spin_bands_to_plot.append(band_data['energy_ev'].tolist())



# Plot bands
for band_data in up_spin_bands_to_plot:
    plt.plot(up_spin_df['|k|'][0:num_k_points], band_data, label=f'Band {band_index + 1}', color='steelblue', linestyle='--', lw=2.0)



# Set custom k-point labels and ticks
plt.xticks(kpoint_ticks, kpoint_labels_plot, fontsize=16)


plt.yticks(fontsize=16)

# plt.xlabel('|k|', fontsize=18)
plt.ylabel('Energy (eV)', fontsize=18)
plt.title(plot_title, fontsize=20)
plt.axhline(y=0, color='black', linestyle='--', linewidth=1.7)  # Add a horizontal line at y=0
# plt.legend()
# Disable horizontal grid lines
plt.grid(which='major', axis='x', linestyle='-', linewidth=1.8, color='gray')
plt.xlim(kpoint_ticks[0], kpoint_ticks[-1])
plt.tight_layout()
plt.show()
