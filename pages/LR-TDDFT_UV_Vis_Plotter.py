import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import interp1d
import io

# Set page config
st.set_page_config(page_title='LR-TDDFT UV/Vis Plotter 📈', layout='wide', page_icon="⚛️",
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


st.title('LR-TDDFT Spectrum Visualizer')
st.write('Upload or paste your Turbomole LR-TDDFT spectrum file')

# Conversion constants
EV_TO_NM = 1239.84
EV_TO_CM = 8065.54
EV_TO_AU = 0.03674932

def detect_input_units(text):
    """Detect units from the input text header"""
    lines = text.strip().split('\n')
    for line in lines:
        line_lower = line.lower()
        if 'cm^(-1)' in line_lower or 'cm-1' in line_lower or '1/cm' in line_lower:
            return '1/cm'
        elif 'hartree' in line_lower or 'a.u.' in line_lower:
            return 'au'
        elif 'nm' in line_lower and 'energy' in line_lower:
            return 'nm'
        elif 'ev' in line_lower:
            return 'eV'
    return 'eV'  # default

def convert_to_ev(energy, from_unit):
    """Convert energy from any unit to eV"""
    if from_unit == 'eV':
        return energy
    elif from_unit == 'nm':
        return EV_TO_NM / energy
    elif from_unit == '1/cm':
        return energy / EV_TO_CM
    elif from_unit == 'au':
        return energy / EV_TO_AU
    return energy

def convert_from_ev(energy_ev, to_unit):
    """Convert energy from eV to any unit"""
    if to_unit == 'eV':
        return energy_ev
    elif to_unit == 'nm':
        return EV_TO_NM / energy_ev
    elif to_unit == '1/cm':
        return energy_ev * EV_TO_CM
    elif to_unit == 'au':
        return energy_ev * EV_TO_AU
    return energy_ev

# File upload or text input
input_method = st.radio('Choose input method:', ['Paste Text', 'Upload File'])

raw_input = None
if input_method == 'Paste Text':
    raw_input = st.text_area('Paste your spectrum data here:', height=300)
else:
    uploaded_file = st.file_uploader('Upload spectrum file', type=['txt', 'dat'])
    if uploaded_file is not None:
        raw_input = uploaded_file.read().decode('utf-8')

if raw_input:
    # Detect input units
    detected_units = detect_input_units(raw_input)
    
    # Allow user to override detected units
    units_col1, units_col2 = st.columns(2)
    input_units = units_col1.selectbox('Input units (auto-detected)', 
                                       ['eV', 'nm', '1/cm', 'au'],
                                       index=['eV', 'nm', '1/cm', 'au'].index(detected_units))
    output_units = units_col2.selectbox('Output units for plotting', 
                                        ['eV', 'nm', '1/cm', 'au'],
                                        index=0)
    
    # Clean the input (remove comment lines)
    lines = raw_input.strip().split('\n')
    cleaned_lines = [line for line in lines if not line.strip().startswith('#')]
    cleaned_input = '\n'.join(cleaned_lines)
    
    # Parse the data
    st.write('### Parsed LR-TDDFT Spectrum')
    df = pd.read_csv(io.StringIO(cleaned_input), delim_whitespace=True, 
                     names=['Energy', 'Intensity'], header=None)
    
    # Show original data
    df_display = df.copy()
    df_display.columns = [f'Energy ({input_units})', 'Intensity']
    st.dataframe(df_display.style.format({f"Energy ({input_units})": "{:.7f}", "Intensity": "{:.7e}"}))
    
    # Convert input energies to eV (internal representation)
    omega0_ev = convert_to_ev(df['Energy'].values, input_units)
    osc_strength0 = df['Intensity'].values
    
    # Convert to output units
    omega0 = convert_from_ev(omega0_ev, output_units)
    
    # Broadening parameters
    st.write('### Broadening Parameters')
    broad_col1, broad_col2, broad_col3 = st.columns(3)
    
    broadening_type = broad_col1.selectbox('Broadening Type', 
                                           ['Gaussian', 'Lorentzian'])
    
    broadening_width_ev = broad_col2.number_input('Broadening Width (eV)', 
                                                  min_value=0.001, max_value=2.0, 
                                                  value=0.15, step=0.01, 
                                                  format="%.3f")
    
    num_points = broad_col3.number_input('Number of Points', 
                                         min_value=100, max_value=10000, 
                                         value=1000, step=100)
    
    # Create broadened spectrum
    minOmega = min(omega0) - abs(0.1 * (max(omega0) - min(omega0)))
    maxOmega = max(omega0) + abs(0.1 * (max(omega0) - min(omega0)))
    
    # Handle reversed scales (nm)
    if minOmega > maxOmega:
        minOmega, maxOmega = maxOmega, minOmega
    
    omega = np.linspace(minOmega, maxOmega, num_points)
    osc_strength = np.zeros_like(omega)
    
    # Convert omega and omega0 to eV for broadening calculation
    omega_ev = convert_to_ev(omega, output_units)
    
    if broadening_type == 'Gaussian':
        for i, (e0_ev, f0) in enumerate(zip(omega0_ev, osc_strength0)):
            osc_strength += f0 * np.exp(-((omega_ev - e0_ev)**2) / (2 * broadening_width_ev**2))
    else:  # Lorentzian
        for i, (e0_ev, f0) in enumerate(zip(omega0_ev, osc_strength0)):
            osc_strength += f0 * (broadening_width_ev / (2 * np.pi)) / \
                           ((omega_ev - e0_ev)**2 + (broadening_width_ev / 2)**2)
    
    # Plotting options
    st.write('### Plotting Options')
    
    # Peak finding
    peaks_col1, peaks_col2 = st.columns(2)
    isFP = peaks_col1.checkbox('Find Peaks', value=False)
    prominence_val = peaks_col2.number_input('Prominence value', 
                                            min_value=0.0, value=0.01, step=0.01)
    
    # Normalization
    normalize_col1, normalize_col2 = st.columns(2)
    isNormalize = normalize_col1.checkbox('Normalize', value=True)
    norm_factor = normalize_col2.number_input('Value by which to normalize the highest peak', 
                                              min_value=0.01, value=1.0, step=0.1)
    
    if isNormalize:
        peak_indices, _ = signal.find_peaks(osc_strength, prominence=prominence_val)
        if len(peak_indices) > 0:
            highest_peak_val = np.max(osc_strength[peak_indices])
            osc_strength = (osc_strength / highest_peak_val) * norm_factor
            osc_strength0 = (osc_strength0 / highest_peak_val) * norm_factor
    
    # Axis limits
    lim_col1, lim_col2, lim_col3, lim_col4 = st.columns(4)
    min_osc_strength = min(osc_strength) - 0.1 * max(osc_strength)
    max_osc_strength = max(osc_strength) + 0.1 * max(osc_strength)
    
    lxlim = lim_col1.number_input('Lower x-axis limit', value=float(minOmega), format="%.4f")
    uxlim = lim_col2.number_input('Upper x-axis limit', value=float(maxOmega), format="%.4f")
    lylim = lim_col3.number_input('Lower y-axis limit', value=float(min_osc_strength), format="%.4f")
    uylim = lim_col4.number_input('Upper y-axis limit', value=float(max_osc_strength), format="%.4f")
    
    # Titles
    title_col1, title_col2, title_col3 = st.columns(3)
    xtitle = title_col1.text_input('X-axis label', value=f'Energy ({output_units})')
    ytitle = title_col2.text_input('Y-axis label', value='Intensity (arb. units)')
    plot_title = title_col3.text_input('Plot title', value='Absorption Spectrum')
    
    # Figure size
    figsize_col1, figsize_col2 = st.columns(2)
    figsize_width = figsize_col1.number_input('Plot width', min_value=5.0, value=10.0, step=0.5)
    figsize_height = figsize_col2.number_input('Plot height', min_value=3.0, value=6.0, step=0.5)
    
    # Spectrum display options
    st.write('### Spectrum Display')
    spectrum_display = st.selectbox('Spectrum Display Mode', 
                                   ['Broadened + Sticks', 'Broadened Only', 'Sticks Only'])
    
    plot_type = st.selectbox('Broadened spectrum style', 
                            ['outline', 'shaded', 'shaded+outline', 'bar'])
    
    # Plot styling
    plot_col1, plot_col2, plot_col3, plot_col4 = st.columns(4)
    plot_color = plot_col1.color_picker('Spectrum color', '#0023F9')
    stick_color = plot_col2.color_picker('Stick spectrum color', '#FF0000')
    transparency = plot_col3.checkbox('Transparent background', value=True)
    line_width = plot_col4.slider('Line Width', 0.4, 10.0, 3.0, step=0.1)
    
    line_style = st.selectbox('Line style', ['solid', 'dashed', 'dotted', 'dashdot'])
    
    # Create plot
    st.write('### Spectrum Plot')
    fig, ax = plt.subplots(figsize=[figsize_width, figsize_height])
    ax.tick_params(axis='both', which='major', labelsize=24)
    ax.tick_params(axis='both', which='minor', labelsize=24)
    
    # Plot based on display mode
    if spectrum_display in ['Broadened + Sticks', 'Broadened Only']:
        if plot_type == 'bar':
            ax.bar(omega, osc_strength, color=plot_color, lw=line_width, width=(omega[1]-omega[0]))
        elif plot_type == 'outline':
            ax.plot(omega, osc_strength, color=plot_color, lw=line_width, linestyle=line_style)
        elif plot_type == 'shaded':
            ax.fill_between(omega, lylim, osc_strength, facecolor=plot_color, alpha=0.5)
        elif plot_type == 'shaded+outline':
            ax.plot(omega, osc_strength, color=plot_color, lw=line_width, linestyle=line_style)
            ax.fill_between(omega, lylim, osc_strength, facecolor=plot_color, alpha=0.5)
    
    if spectrum_display in ['Broadened + Sticks', 'Sticks Only']:
        markerline, stemlines, baseline = ax.stem(omega0, osc_strength0, linefmt=stick_color, 
                                                   markerfmt='o', basefmt=' ')
        plt.setp(stemlines, 'linewidth', 1.5)
        plt.setp(markerline, 'markersize', 0)
    
    ax.set_xlabel(xtitle, fontsize=24)
    ax.set_ylabel(ytitle, fontsize=24)
    ax.set_xlim([lxlim, uxlim])
    ax.set_ylim([lylim, uylim])
    ax.set_title(plot_title, fontsize=27)
    
    if isFP:
        peak_indices, _ = signal.find_peaks(osc_strength, prominence=prominence_val)
        for i in peak_indices:
            ax.annotate(f'({omega[i]:.4f}, {osc_strength[i]:.4f})', 
                       [omega[i], osc_strength[i]], fontsize=10)
    
    plt.tight_layout()
    plt.savefig('spectrum.png', transparent=transparency, dpi=300)
    st.pyplot(fig)
    
    with open('spectrum.png', 'rb') as f:
        st.download_button('Download plot as PNG', f, file_name='spectrum.png')
    
    # Color calculation
    st.write('### Emission/Absorption Color')
    
    def wavelength_to_rgb(wavelength):
        """Convert wavelength (nm) to RGB"""
        if wavelength < 380 or wavelength > 750:
            return (0, 0, 0)
        
        if 380 <= wavelength < 440:
            R = -(wavelength - 440) / (440 - 380)
            G = 0.0
            B = 1.0
        elif 440 <= wavelength < 490:
            R = 0.0
            G = (wavelength - 440) / (490 - 440)
            B = 1.0
        elif 490 <= wavelength < 510:
            R = 0.0
            G = 1.0
            B = -(wavelength - 510) / (510 - 490)
        elif 510 <= wavelength < 580:
            R = (wavelength - 510) / (580 - 510)
            G = 1.0
            B = 0.0
        elif 580 <= wavelength < 645:
            R = 1.0
            G = -(wavelength - 645) / (645 - 580)
            B = 0.0
        else:
            R = 1.0
            G = 0.0
            B = 0.0
        
        return (int(R * 255), int(G * 255), int(B * 255))
    
    # Calculate weighted average wavelength (convert to nm for color calculation)
    omega_nm = convert_from_ev(omega0_ev, 'nm')
    
    valid_visible = (omega_nm >= 380) & (omega_nm <= 750)
    
    if np.any(valid_visible):
        weighted_wl = np.average(omega_nm[valid_visible], 
                                weights=osc_strength0[valid_visible])
        color_rgb = wavelength_to_rgb(weighted_wl)
        color_hex = '#{:02x}{:02x}{:02x}'.format(*color_rgb)
        
        col1, col2 = st.columns(2)
        col1.write(f'**Dominant wavelength:** {weighted_wl:.1f} nm')
        col2.markdown(f'**Color:** <span style="background-color:{color_hex};padding:10px 30px;">&nbsp;&nbsp;&nbsp;&nbsp;</span>', 
                     unsafe_allow_html=True)
    else:
        st.write('No absorption in visible range (380-750 nm)')
    
    # Download broadened spectrum
    st.write('### Download Broadened Spectrum')
    chart_data = pd.DataFrame({
        f'Energy ({output_units})': omega,
        'Intensity': osc_strength
    })
    
    st.dataframe(chart_data)
    
    export_spec_file_name = 'lr-tddft_spectrum_broadened.xlsx'
    chart_data.to_excel(export_spec_file_name, index=False)
    with open(export_spec_file_name, 'rb') as f:
        st.download_button('Download as EXCEL', f, file_name=export_spec_file_name)

else:
    st.info('Please paste or upload your spectrum data to begin.')
