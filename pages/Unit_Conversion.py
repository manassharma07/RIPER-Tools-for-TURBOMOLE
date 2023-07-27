import streamlit as st
import pandas as pd

# Sidebar stuff
st.sidebar.write('# About')
st.sidebar.write('### Originally Made By [Manas Sharma](https://manas.bragitoff.com)')
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

st.write('## Unit Conversion Tool')
st.write('This online tool allows you to inter-convert between various units mainly used in RT-TDDFT calculations via RIPER.')

# Supported categories
categories = ['energy', 'length', 'time', 'frequency', 'elec']
# Supported units for energy
energyUnits = ['ev', 'au', 'ha', 'ryd', 'kj/mol', 'kcal/mol', 'kcal', 'joule', 'kj']
# Supported units for distance
distUnits = ['bohr', 'au', 'angs', 'um', 'nm', 'pm', 'fm']
# Supported units for time
timeUnits = ['au', 'fs', 'ps', 'as', 's']
# Supported units for frequency/wavelength/time period
freqUnits = ['au', 'ev', 'cm^-1', 'hz', 'mhz', 'thz', 'nm', 'um']
# Supported units for electric field amplitude
elecUnits = ['au', 'v/cm', 'v/m']

# Dictionary of supported units and categories
catUnitDict = {'energy': energyUnits, 'length': distUnits, 'frequency': freqUnits, 'time': timeUnits, 'elec': elecUnits}

# Conversion factor from arbitrary unit to au for Energy
energy2au = {'ev': 0.036749405469679, 'au': 1.0, 'ha': 1.0, 'ryd': 0.5, 'kj/mol': 3.8088E-4, 'kcal/mol': 1.5936E-3,
            'kcal': 9.5968937625695E+20, 'joule': 2.2937126583579E+17, 'kj': 2.2937126583579E+20}
# Conversion factor from arbitrary unit to au for Distance
dist2au = {'bohr': 1.0, 'au': 1.0, 'angs': 1.8897268777744, 'um': 18897.268777744, 'nm': 18.897268777744,
        'pm': 0.018897268777744, 'fm': 1.8897268777744E-5}
# Conversion factor from arbitrary unit to au for Time
time2au = {'au': 1.0, 'fs': 41.341374575751, 'ps': 41341.374575751, 'as': 0.041341374575751, 's': 4.1341374575751E+16}
# Conversion factor from arbitrary unit to au for Frequency
freq2au = {'au': 1.0, 'ev': 0.036749405469679, 'cm^-1': 4.556335270832134e-06, 'hz': 1.519828500716E-16,
        'mhz': 1.519828500716E-10, 'thz': 1.519828500716E-4}
# Conversion factor from arbitrary unit to au for Electric field amplitude
elec2au = {'au': 1.0, 'v/cm': 0.19446898091988207E-9, 'v/m': 0.19446898091988207E-11}
# Conversion factor to go from frequency in au to wavelength in
freq2wlau = {'nm': 45.563453803879376, 'um': 0.045563453803879376}
# Conversion factor to go from frequency in au to time period in
freq2tpau = {'fs': 1.519828500716E-1, 'ps': 1.519828500716E-4}
# Conversion factor to go from electric field amplitude in au to intensity in
elec2inau = {'w/cm^2': 3.509447584910074703E16, 'tw/cm^2': 3.509447584910074703E4}
# Conversion factor to go from time period in au to frequency in
tp2freqau = {'ev': 170.97364851478514}
# Conversion factor to go from time period in au to wavelength in
tp2wlau = {'nm': 7.251655098725834}
# Conversion factor to go from wavelength in au to frequency in
wl2freqau = {'ev': 23429.626193126864}
# Conversion factor to go from wavelength in au to time period in
wl2tpau = {'fs': 0.00017651389861382622}

# Dictionary of conversion factors and categories
catCFDict = {'energy': energy2au, 'length': dist2au, 'frequency': freq2au, 'time': time2au, 'elec': elec2au}


input_type = st.selectbox('Category of units', categories)
input_unit = st.selectbox('Unit of the input value', catUnitDict[input_type])
input_value = st.number_input('Input value to convert', value=1.0, step=None, format="%0.14f")

# Convert the input value to the au of the corresponding category
au_value = catCFDict[input_type][input_unit] * input_value

# Print all the possible conversions
# Create a DataFrame to store the conversions
conversions = []
for unit, factor in catCFDict[input_type].items():
    converted_value = au_value / factor
    conversions.append([converted_value, unit])

# Additional conversions based on input type
if input_type == 'frequency':
    for unit, factor in freq2wlau.items():
        converted_value = factor / au_value
        conversions.append([converted_value, unit])
    for unit, factor in freq2tpau.items():
        converted_value = factor / au_value
        conversions.append([converted_value, unit])

if input_type == 'time':
    for unit, factor in tp2freqau.items():
        converted_value = factor / au_value
        conversions.append([converted_value, unit])
    for unit, factor in tp2wlau.items():
        converted_value = au_value * factor
        conversions.append([converted_value, unit])

if input_type == 'length':
    for unit, factor in wl2freqau.items():
        converted_value = factor / au_value
        conversions.append([converted_value, unit])
    for unit, factor in wl2tpau.items():
        converted_value = au_value * factor
        conversions.append([converted_value, unit])

if input_type == 'elec':
    for unit, factor in elec2inau.items():
        converted_value = factor * (au_value ** 2)
        conversions.append([converted_value, unit])

# Create a DataFrame from the conversions list
df = pd.DataFrame(conversions, columns=['Converted Value', 'Unit'])


# Format large numbers without exponents in the DataFrame
df['Converted Value'] = df['Converted Value'].apply('{:.15g}'.format)

# Display the DataFrame as a table
st.subheader('Conversions:')
st.dataframe(df.set_index('Unit'), width=400)
