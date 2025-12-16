import streamlit as st
import pandas as pd
from scipy import signal
import numpy as np
import io
import matplotlib.pyplot as plt

# Set page config
st.set_page_config(page_title='RT-TDDFT Absorption Spectrum Plotter 📈', layout='wide', page_icon="⚛️",
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

def parse_rtspectrum_content(content):
    lines = content.split('\n')\
    # Remove empty lines at the beginning and end
    lines = [line.strip() for line in lines if line.strip()]

    energy_units = lines[1].split('(')[1].split(')')[0].strip()
    data_lines = [line.split() for line in lines[2:-1]]
    data = pd.DataFrame(data_lines, columns=['Energy (' + energy_units + ')', 'Intensity'])
    data['Energy (' + energy_units + ')'] = data['Energy (' + energy_units + ')'].astype(float)
    data['Intensity'] = data['Intensity'].astype(float)
    
    return data, energy_units


st.title('RT-TDDFT Absorption Spectrum Plotter')
st.write('This utility helps you plot the RT-TDDFT spectrum as saved in the `rtspec` file after running an Absorption Spectrum calculation using `RIPER`.')

#st.sidebar.write('### Input Data')
uploaded_file = st.file_uploader('Upload the `rtspec` file')
pasted_content = st.text_area('or paste the contents of the `rtspec` file here:', height=400)

content = None
osc_strength = []
omega = []

if uploaded_file is not None:
    content = uploaded_file.read().decode('utf-8')
elif pasted_content:
    content = pasted_content

if content:
    st.write('### Parsed data')
    data, energy_units = parse_rtspectrum_content(content)

    st.dataframe(data)

    omega = data['Energy (' + data.columns[0].split('(')[1].split(')')[0].strip() + ')'].values
    osc_strength = data['Intensity'].values


    ## Plot the absorption spectrum
    st.write('### Plotting the RT-TDDFT Spectrum')

    peaks_col1, peaks_col2 = st.columns(2)
    isFP = peaks_col1.checkbox('Find Peaks', value=False)
    prominence_val = peaks_col2.text_input(label='Prominence value', value='0', key='prominence')
    prominence_val = float(prominence_val)

    # omega= 1239.84193 / omega
    # st.write(omega)
    if isFP:
        peak_indices, dict_peak = signal.find_peaks(osc_strength, prominence=prominence_val)
        highest_peak_val = np.max(osc_strength[peak_indices])

    normalize_col1, normalize_col2 = st.columns(2)
    isNormalize = normalize_col1.checkbox('Normalize', value=True)
    norm_factor = normalize_col2.text_input(label='Value by which to normalize the highest peak', value='1', key='norm_factor')
    norm_factor = float(norm_factor)

    if isNormalize:
        peak_indices, dict_peak = signal.find_peaks(osc_strength, prominence=prominence_val)
        highest_peak_val = np.max(osc_strength[peak_indices])
        # st.write(highest_peak_val)
        osc_strength = (osc_strength/highest_peak_val)*norm_factor
        # st.write(np.max(osc_strength))


    lim_col1, lim_col2, lim_col3, lim_col4 = st.columns(4)
    lxlim = lim_col1.text_input(label='Enter the lower limit for x-axis', value=str(min(omega)),key='lxlim')
    lxlim = float(lxlim)
    uxlim = lim_col2.text_input(label='Enter the upper limit for x-axis', value=str(max(omega)),key='uxlim')
    uxlim = float(uxlim)
    # --- Mask data inside x-limits ---
    xmask = (omega >= lxlim) & (omega <= uxlim)

    # Fallback in case the range is empty
    if np.any(xmask):
        osc_in_range = osc_strength[xmask]
        ymin = np.min(osc_in_range)
        ymax = np.max(osc_in_range)
    else:
        ymin = np.min(osc_strength)
        ymax = np.max(osc_strength)

    # Add padding (10%)
    ypad = 0.1 * (ymax - ymin if ymax > ymin else ymax)
    lylim = lim_col3.text_input(label='Enter the lower limit for y-axis', value=str(ymin - ypad))
    lylim = float(lylim)
    uylim = lim_col4.text_input(label='Enter the upper limit for y-axis', value=str(ymax + ypad))
    uylim = float(uylim)

    title_col1, title_col2, title_col3 = st.columns(3)
    xtitle = title_col1.text_input(label='Enter the label for x-axis', value='Energy ('+energy_units+')', key='xtitle')
    ytitle = title_col2.text_input(label='Enter the label for y-axis', value='Intensity (arb. units)', key='ytitle')
    plot_title = title_col3.text_input(label='Enter the title for the plot', value='Absorption Spectrum', key='plot_title')

    figsize_col1, figsize_col2 = st.columns(2)
    figsize_width = figsize_col1.text_input(label='Enter the plot figure width', value='10', key='figsize_width')
    figsize_width = float(figsize_width)
    figsize_height = figsize_col2.text_input(label='Enter the plot figure height', value='6', key='figsize_height')
    figsize_height = float(figsize_height)

    plot_type = st.selectbox('Select the plot style',
        ( 'outline','shaded', 'shaded+outline'))

    plot_col1, plot_col2, plot_col3, plot_col4 = st.columns(4)
    plot_color = plot_col1.color_picker('Pick a Color for the Spectrum plot', '#0023F9')
    transparency = plot_col2.checkbox('Make the plot transparent', value=True)
    line_width = plot_col3.slider('Line Width', 0.4, 10., 3.0, step=0.1)
    line_style = plot_col4.selectbox('Select the line style',
        ( 'solid','dashed', 'dotted', 'dashdot'))


    import matplotlib.pylab as pylab
    # params = {'legend.fontsize': 'xx-large',
    #      'figure.figsize': (figsize_width, figsize_height),
    #     'axes.labelsize': 'xx-large',
    #     'axes.titlesize':'xx-large',
    #     'xtick.labelsize':'xx-large',
    #     'ytick.labelsize':'xx-large'}
    # pylab.rcParams.update(params)

    fig, ax = plt.subplots(figsize=[figsize_width, figsize_height])
    # We change the fontsize of minor ticks label 
    ax.tick_params(axis='both', which='major', labelsize=24)
    ax.tick_params(axis='both', which='minor', labelsize=24)
    # ax.plot(omega, osc_strength, color=plot_color, lw=line_width, linestyle=line_style)
    if plot_type=='outline':
        ax.plot(omega, osc_strength, color=plot_color, lw=line_width, linestyle=line_style)
    if plot_type=='shaded':
        ax.fill_between(omega, lylim, osc_strength, facecolor=plot_color, alpha=0.5)
    if plot_type=='shaded+outline':
        ax.plot(omega, osc_strength, color=plot_color, lw=line_width, linestyle=line_style)
        ax.fill_between(omega, lylim, osc_strength, facecolor=plot_color, alpha=0.5)
    ax.set_title('Absorption Spectrum', fontsize=50)
    ax.set_xlabel(xtitle, fontsize=24)
    ax.set_ylabel(ytitle, fontsize=24)
    ax.set_xlim([lxlim, uxlim])
    ax.set_ylim([lylim, uylim])

    ax.set_title(plot_title, fontsize=27)
    if isFP:
        for i in peak_indices:
            plt.annotate('('+ str(np.round(omega[i],4))+', '+str(np.round(osc_strength[i],4))+' )', [omega[i], osc_strength[i]])

    plt.tight_layout()
    plt.savefig('spectrum.png', transparent=transparency)
    st.pyplot(fig)
    with open('spectrum.png', 'rb') as f:
        st.download_button('Download plot as PNG file', f, file_name='spectrum.png')  

   
