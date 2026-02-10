from fileinput import filename
import streamlit as st
import io  
import pandas as pd 
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
from scipy.fft import fft
import os
from numba import njit, jit

# Set page config
st.set_page_config(page_title='HHG Spectrum Calculator and Plotter', layout='wide', page_icon="⚛️",
menu_items={
         'About': "A web app to help you with DFT related calculations using the RIPER module of [TURBOMOLE](https://www.turbomole.org/)"
     })

# Sidebar stuff
st.sidebar.write('# About')
st.sidebar.write(' Made By [Manas Sharma](https://manas.bragitoff.com)')
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


st.title('Turbomole RT-TDDFT High Harmonic Generation (HHG) Spectrum Calculator')
st.write('This tool allows you to calculate the HHG spectra by entering the contents of the `rtdipo` file as well as the `control` file.')


# read the sample rtdipo_HHG_Sample file and put its contents in a string variable to be used as a placeholder in the text area for rtdipo input
placeholder_rtdipo_str = ''
with open('pages/rtdipo_HHG_Sample', 'r') as f:
     placeholder_rtdipo_str = f.read()

# Supported units for energy
energy_units = ['ev','au','ha','ryd','kj/mol','kcal/mol','kcal','joule','kj']
# Conversion factor from arbitrary unit to au for Energy
energy2au = {'ev':0.036749405469679, 'au':1.0, 'ha':1.0, 'ryd':0.5, 'kj/mol':3.8088E-4, 'kcal/mol':1.5936E-3, \
'kcal':9.5968937625695E+20, 'joule':2.2937126583579E+17, 'kj':2.2937126583579E+20 }


time_step = st.text_input(label = 'Enter the time step (a.u.)', value='0.5')
time_step = float(time_step)

## Electric field details
st.write('### Electric Field Details')
ex_col1, ex_col2, ex_col3 = st.columns(3)
amp_x = ex_col1.text_input(label = 'Enter the amplitude along x (a.u.)', value='0.0',key='Ex')
amp_x = float(amp_x)
amp_y = ex_col2.text_input(label = 'Enter the amplitude along y (a.u.)', value='0.00653769',key='Ey')
amp_y = float(amp_y)
amp_z = ex_col3.text_input(label = 'Enter the amplitude along z (a.u.)', value='0.0',key='Ez')
amp_z = float(amp_z)
phase_col1, phase_col2, phase_col3 = st.columns(3)
phase_x = phase_col1.text_input(label = 'Enter the phase along x (radians)', value='0.0',key='phase_x')
phase_x = float(phase_x)
phase_y = phase_col2.text_input(label = 'Enter the phase along y (radians)', value='0.0',key='phase_y')
phase_y = float(phase_y)
phase_z = phase_col3.text_input(label = 'Enter the phase along z (radians)', value='0.0',key='phase_z')
phase_z = float(phase_z)
E_par_col1, E_par_col2 = st.columns(2)
omega0_E = E_par_col1.text_input(label = 'Enter the fundamental frequency of the laser pulse (a.u.)', value='0.0227816',key='omega0_E')
omega0_E = float(omega0_E)
sigma_E = E_par_col2.text_input(label = 'Enter the FWHM of the laser pulse (a.u.)', value='1241.099484208',key='sigma_fwhm_E')
sigma_E = float(sigma_E)
# Total duration of the pulse
pulse_duration = 2*sigma_E



st.write('### INPUT (rtdipo)')
input_rtdipo_str = st.text_area(label='Enter the contents of `rtdipo` here', value = placeholder_rtdipo_str, placeholder = 'Put your text here', height=400)

def clean_input_str(input_rtdipo_str):
     cleaned_lines = []
     for line in input_rtdipo_str.splitlines():
          # Remove comments (# and $)
          line = line.partition('#')[0].partition('$')[0].rstrip()
          # Keep non-empty lines
          if line:
               cleaned_lines.append(line)
     return os.linesep.join(cleaned_lines)



cleaned_input = clean_input_str(input_rtdipo_str)




# pd.set_option("display.precision", 14)
st.write('### Parsed dipole moments along x, y and z at different timesteps')
# df = pd.read_csv(io.StringIO(cleaned_input), delim_whitespace=True, names=['Time step', 'x direction', 'y direction', 'z direction'])
# df = pd.read_csv(io.StringIO(cleaned_input), index_col = False,sep="   ", header=0,names=['Time step', 'x direction', 'y direction', 'z direction'])
df = pd.read_csv(io.StringIO(cleaned_input), skiprows=1, delim_whitespace=True, names=['Time step', 'x direction', 'y direction', 'z direction'])
# Insert the actual time column in the dataframe object
df.insert(1,'Time (a.u.)', time_step*df['Time step'])
# df = df.set_index('Time step')
# st.dataframe(df.style.format("{:.8%}"))
# Set precision for printing
st.dataframe(df.style.format({"E": "{:.2f}"}))
# st.dataframe(df)
# st.write(input_rtdipo_str)


# determining the name of the file
export_rtdipo_file_name = 'rtdipo.xlsx'

# saving to excel file
df.to_excel(export_rtdipo_file_name)
with open(export_rtdipo_file_name, 'rb') as f:
     st.download_button('Download EXCEL file', f, file_name=export_rtdipo_file_name)  
# st.write('rtdipo is written to Excel File successfully.')


# Plotting the dipole moments
st.write('### Plots of the dipole moments')
dip_col1, dip_col2, dip_col3 = st.columns(3)
dip_col1.write('Dipole moment along x: $\mu_x$')
fig, ax = plt.subplots()
ax.plot(df['Time (a.u.)'], df['x direction'])
ax.set_title('Dipole moment along x-axis', fontsize=16)
ax.set_xlabel('Time (a.u.)', fontsize=14)
ax.set_ylabel('Dipole Moment (a.u.)', fontsize=14)
dip_col1.pyplot(fig)

dip_col2.write('Dipole moment along y: $\mu_y$')
fig, ax = plt.subplots()
ax.plot(df['Time (a.u.)'], df['y direction'])
ax.set_title('Dipole moment along y-axis', fontsize=16)
ax.set_xlabel('Time (a.u.)', fontsize=14)
ax.set_ylabel('Dipole Moment (a.u.)', fontsize=14)
dip_col2.pyplot(fig)

dip_col3.write('Dipole moment along z: $\mu_z$')
fig, ax = plt.subplots()
ax.plot(df['Time (a.u.)'], df['z direction'])
ax.set_title('Dipole moment along z-axis', fontsize=16)
ax.set_xlabel('Time (a.u.)', fontsize=14)
ax.set_ylabel('Dipole Moment (a.u.)', fontsize=14)
dip_col3.pyplot(fig)

## Start the spectra calculation
# Create a time vector
time = df['Time (a.u.)']
time = np.array(time)

# Create a dipole moment matrix
mu = np.zeros((len(df['Time step']), 3))
mu[:,0] = df['x direction']
mu[:,1] = df['y direction']
mu[:,2] = df['z direction']

# Create an electric field matrix 
@njit
def f_envelope(time, sigma):
   value = np.zeros(time.shape[0])
   for i in range(time.shape[0]):
       t = time[i]
       if t>0 and t<2*sigma:
           value[i]=np.cos(np.pi/2/sigma*(t-sigma))**2
   return value

E = np.zeros((len(df['Time step']), 3))
E[:,0] = amp_x*np.sin(omega0_E*time+phase_x)*f_envelope(time, sigma_E)
E[:,1] = amp_y*np.sin(omega0_E*time+phase_y)*f_envelope(time, sigma_E)
E[:,2] = amp_z*np.sin(omega0_E*time+phase_z)*f_envelope(time, sigma_E)
ampE = np.array([amp_x, amp_y, amp_z])


# Plot the electric field 
st.write('### Plots of the electric field')
E_col1, E_col2, E_col3 = st.columns(3)
E_col1.write('Electric field along x: $E_x$')
fig, ax = plt.subplots()
ax.plot(time, E[:,0])
ax.set_title('Electric field along x-axis', fontsize=16)
ax.set_xlabel('Time (a.u.)', fontsize=14)
ax.set_ylabel('Electric field (a.u.)', fontsize=14)
E_col1.pyplot(fig)

E_col2.write('Electric field along y: $E_y$')
fig, ax = plt.subplots()
ax.plot(time, E[:,1])
ax.set_title('Electric field along y-axis', fontsize=16)
ax.set_xlabel('Time (a.u.)', fontsize=14)
ax.set_ylabel('Electric field (a.u.)', fontsize=14)
E_col2.pyplot(fig)

E_col3.write('Electric field along z: $E_z$')
fig, ax = plt.subplots()
ax.plot(time, E[:,2])
ax.set_title('Electric field along z-axis', fontsize=16)
ax.set_xlabel('Time (a.u.)', fontsize=14)
ax.set_ylabel('Electric field (a.u.)', fontsize=14)
E_col3.pyplot(fig)

## Frequency details for transformation
st.write('### Frequency/Energy Details for RT-TDDFT Spectrum')
w_col1, w_col2, w_col3 = st.columns(3)
omega_min = w_col1.text_input(label = 'Enter the lower limit of frequency/energy (a.u.)', value='0.0',key='omega_min')
omega_min = float(omega_min)
omega_max = w_col2.text_input(label = 'Enter the upper limit of frequency/energy (a.u.)', value='1.2',key='omega_max')
omega_max = float(omega_max)
d_omega = w_col3.text_input(label = 'Enter the frequency/energy spacing (a.u.)', value='0.001',key='d_omega')
d_omega = float(d_omega)

energy_units = st.selectbox('Select the units for energy',
     energy_units)
unit_conv_factor = 1./energy2au[energy_units]
c = 137.036

## Fourier Transformation details
st.write('### Fourier Transformation details for RT-TDDFT HHG Spectrum')
# Max. time that should be used for Fourier Transformation (spectrum calculation)
max_time = st.text_input(label = 'Enter the max. value of time that should be used for Fourier transform (Cannot be more than the propagation time)', value=str(time[-1]),key='max_time') # a.u.
max_time = float(max_time)
max_time_step = int(max_time/time_step)-1
omega = np.arange(omega_min, omega_max, d_omega)

if not np.isclose(amp_x,0.0):
   x_comp = mu[0:max_time_step,0]
else:
   x_comp = np.zeros(max_time_step)
if not np.isclose(amp_y,0.0):
   y_comp = mu[0:max_time_step,1]
else:
   y_comp = np.zeros(max_time_step)
if not np.isclose(amp_z,0.0):
   z_comp = mu[0:max_time_step,2]
else:
   z_comp = np.zeros(max_time_step)

isSavgol = st.checkbox('Savgol', value=False)
gauge = st.selectbox('Select the gauge',
     ( 'length','velocity', 'acceleration' ))
if gauge=='length':
    isAG=False
    isVG=False
if gauge=='velocity':
    isVG=True
    isAG=False
if gauge=='acceleration':
    isVG=False
    isAG=True

# Use smoothing/damping for smooth termination of dipole moment towards the end of the pulse
isSmoothingOn = st.checkbox('Smoothing', value=False)
if isSmoothingOn:
   smoothing_type =  st.selectbox('Select the smoothing type',
     ( 'Cossq','Sinsq', 'Falvo' ))
   if smoothing_type=='Sinsq':
       isSinsqSmoothing = True 
       isCossqSmoothing = False
       isFalvosqSmoothing = False
   if smoothing_type=='Cossq':
       isSinsqSmoothing = False 
       isCossqSmoothing = True
       isFalvosqSmoothing = False
   if smoothing_type=='Falvo':
       isSinsqSmoothing = False 
       isCossqSmoothing = False
       isFalvosqSmoothing = True
   # The tau parameter defines the time value till which 
   # no smoothing/damping should be applied. Ex: tau = 0.9*T
   # means that the smoothing will only be applied to the 
   # dipole moment during the last 10% of the laser pulse.
   tau = st.text_input(label = 'Enter the value of tau such that the smoothing is applied to (1-tau)*Total_pulse_time', value='0.75',key='tau')
   tau = float(d_omega)*pulse_duration
   # The time step (int) from which smoothing/damping is required
   tau_step = int(tau/time_step)
# Use window functions for removing spectral leakage from Fourier transform
isWindowing = st.checkbox('Use windowing', value=True)
windowType = 'hann'

# Use moving averages
isMovingAverages = st.checkbox('Use moving averages', value=False)
if isMovingAverages:
   windowSize = st.text_input(label = 'Enter the window size for moving averages (a.u.)', value='0.003',key='windowSize')
   windowSize = float(windowSize)



def moving_average (values, window_size):
   weights = np.ones(window_size)/window_size
   sma = np.convolve(values, weights, 'same')
   return sma


if isSavgol:
   x_comp = signal.savgol_filter(x_comp, 33, 2, mode='nearest')
   y_comp = signal.savgol_filter(y_comp, 9, 2)
   z_comp = signal.savgol_filter(z_comp, 33, 2)

if isVG:
   x_comp = np.gradient(x_comp, time_step) # First derivative
   y_comp = np.gradient(y_comp, time_step) # First derivative
   z_comp = np.gradient(z_comp, time_step) # First derivative

if isAG:
   x_comp = np.gradient(x_comp, time_step) # First derivative
   y_comp = np.gradient(y_comp, time_step) # First derivative
   z_comp = np.gradient(z_comp, time_step) # First derivative

   x_comp = np.gradient(x_comp, time_step) # Second derivative
   y_comp = np.gradient(y_comp, time_step) # Second derivative
   z_comp = np.gradient(z_comp, time_step) # Second derivative

P = np.zeros(omega.shape[0]) # HHG Power spectrum as a function of omega
Intensity = np.zeros(omega.shape[0]) # log10(P) (a function of omega)
ti = time[0:max_time_step]       

# Smooth termination of dipole moment
if isSmoothingOn:
   # Smoothing is only applied to the dipole moment during the last T-tau duration of the pulse
   if isSinsqSmoothing:
       x_comp[tau_step:] = x_comp[tau_step:]*(np.sin(np.pi*(pulse_duration-ti[tau_step:])/(2*pulse_duration)))**2
       y_comp[tau_step:] = y_comp[tau_step:]*(np.sin(np.pi*(pulse_duration-ti[tau_step:])/(2*pulse_duration)))**2
       z_comp[tau_step:] = z_comp[tau_step:]*(np.sin(np.pi*(pulse_duration-ti[tau_step:])/(2*pulse_duration)))**2
   if isCossqSmoothing:
       # This is actually same as before (sinsq)
       x_comp[tau_step:] = x_comp[tau_step:]*np.cos( np.pi/2 * (ti[tau_step:]-tau)/(pulse_duration-tau) )**2
       y_comp[tau_step:] = y_comp[tau_step:]*np.cos( np.pi/2 * (ti[tau_step:]-tau)/(pulse_duration-tau) )**2
       z_comp[tau_step:] = z_comp[tau_step:]*np.cos( np.pi/2 * (ti[tau_step:]-tau)/(pulse_duration-tau) )**2
   if isFalvosqSmoothing:
       x_comp[tau_step:] = x_comp[tau_step:]*np.exp(-0.5 * ((ti[tau_step:]-tau)/((pulse_duration-tau)/3.0))**4)
       y_comp[tau_step:] = y_comp[tau_step:]*np.exp(-0.5 * ((ti[tau_step:]-tau)/((pulse_duration-tau)/3.0))**4)
       z_comp[tau_step:] = z_comp[tau_step:]*np.exp(-0.5 * ((ti[tau_step:]-tau)/((pulse_duration-tau)/3.0))**4)

#Hann window function
if isWindowing:
   windows = signal.get_window(windowType, max_time_step)

   x_comp = x_comp*windows
   y_comp = y_comp*windows
   z_comp = z_comp*windows


i=0
for omegai in omega:
   # sumP = 0.0
   sumP1 = np.array([0.0, 0.0, 0.0])
   sumP2 = np.array([0.0, 0.0, 0.0])
   temp1 = np.sin(omegai*ti)#*time_step
   temp2 = np.cos(omegai*ti)#*time_step
   sumP1[0] = np.trapz(x_comp*temp1, dx=ti[1]-ti[0]) 
   sumP2[0] = np.trapz(x_comp*temp2, dx=ti[1]-ti[0]) 
   sumP1[1] = np.trapz(y_comp*temp1, dx=ti[1]-ti[0]) 
   sumP2[1] = np.trapz(y_comp*temp2, dx=ti[1]-ti[0]) 
   sumP1[2] = np.trapz(z_comp*temp1, dx=ti[1]-ti[0]) 
   sumP2[2] = np.trapz(z_comp*temp2, dx=ti[1]-ti[0]) 
   P[i]=sumP1[0]**2 + sumP2[0]**2 + sumP1[1]**2 + sumP2[1]**2 + sumP1[2]**2 + sumP2[2]**2
   i=i+1


if isMovingAverages:
   window_size = int(windowSize/ d_omega)
   P = moving_average(P, window_size)
Intensity = np.log10(P)



osc_strength = Intensity
omega = omega*unit_conv_factor






## Plot the absorption spectrum
st.write('### Plotting the RT-TDDFT Spectrum')

peaks_col1, peaks_col2 = st.columns(2)
isFP = peaks_col1.checkbox('Find Peaks', value=False)
prominence_val = peaks_col2.text_input(label='Prominence value', value='0', key='prominence')
prominence_val = float(prominence_val)


if isFP:
     peak_indices, dict_peak = signal.find_peaks(osc_strength, prominence=prominence_val)
     highest_peak_val = np.max(osc_strength[peak_indices])

normalize_col1, normalize_col2 = st.columns(2)
isNormalize = normalize_col1.checkbox('Normalize', value=True)
norm_factor = normalize_col2.text_input(label='Value by which to normalize the highest peak', value='0', key='norm_factor')
norm_factor = float(norm_factor)

if isNormalize:
     peak_indices, dict_peak = signal.find_peaks(osc_strength, prominence=prominence_val)
     highest_peak_val = np.max(osc_strength[peak_indices])
     # st.write(highest_peak_val)
     osc_strength = osc_strength/highest_peak_val
     # st.write(np.max(osc_strength))


lim_col1, lim_col2, lim_col3, lim_col4 = st.columns(4)
lxlim = lim_col1.text_input(label='Enter the lower limit for x-axis', value=str(min(omega)),key='lxlim')
lxlim = float(lxlim)
uxlim = lim_col2.text_input(label='Enter the upper limit for x-axis', value=str(max(omega)),key='uxlim')
uxlim = float(uxlim)
lylim = lim_col3.text_input(label='Enter the lower limit for y-axis', value=str(min(osc_strength) - 0.1*max(osc_strength)),key='lylim')
lylim = float(lylim)
uylim = lim_col4.text_input(label='Enter the upper limit for y-axis', value=str(max(osc_strength) + 0.1*max(osc_strength)),key='uylim')
uylim = float(uylim)

title_col1, title_col2, title_col3 = st.columns(3)
xtitle = title_col1.text_input(label='Enter the label for x-axis', value='Energy ('+energy_units+')', key='xtitle')
ytitle = title_col2.text_input(label='Enter the label for y-axis', value='Intensity (arb. units) \n $\log_{10}|P(\omega)|$', key='ytitle')
plot_title = title_col3.text_input(label='Enter the title for the plot', value="HHG Spectrum (Oscillations along y-axis) \n  time 2758 a.u.", key='plot_title')

figsize_col1, figsize_col2 = st.columns(2)
figsize_width = figsize_col1.text_input(label='Enter the plot figure width', value='10', key='figsize_width')
figsize_width = float(figsize_width)
figsize_height = figsize_col2.text_input(label='Enter the plot figure height', value='6', key='figsize_height')
figsize_height = float(figsize_height)

plot_col1, plot_col2, plot_col3, plot_col4 = st.columns(4)
plot_color = plot_col1.color_picker('Pick a Color for the Spectrum plot', '#0023F9')
transparency = plot_col2.checkbox('Make the plot transparent', value=True)
line_width = plot_col3.slider('Line Width', 0.4, 10., 1.5, step=0.1)
isGrids = plot_col4.checkbox('Grids', value=False)





fig, ax = plt.subplots(figsize=[figsize_width, figsize_height])
 # We change the fontsize of minor ticks label 
ax.tick_params(axis='both', which='major', labelsize=22)
ax.tick_params(axis='both', which='minor', labelsize=22)
ax.plot(omega, osc_strength, color=plot_color, lw=line_width)
ax.set_title('Absorption Spectrum')
ax.set_xlabel(xtitle, fontsize=22)
ax.set_ylabel(ytitle, fontsize=22)
ax.set_xlim([lxlim, uxlim])
ax.set_ylim([lylim, uylim])
ax.set_title(plot_title, fontsize=27)
ax2 = ax.twiny()   # mirror the current x axis
# ax2.tick_params(axis='both', which='minor', labelsize=22)
if isGrids:
   ax2.grid(True)
ax2.set_xlabel('Harmonic Order $\omega/\omega_0$', fontsize=22)
ax2.set_xticks(np.arange(float(lxlim), float(uxlim), step=omega0_E*unit_conv_factor),np.arange(0,int(float(uxlim)/(omega0_E*unit_conv_factor))+1,1), fontsize=16)
ax2.set_xlim([lxlim, uxlim])


if isFP:
     for i in peak_indices:
          plt.annotate('('+ str(np.round(omega[i],4))+', '+str(np.round(osc_strength[i],4))+' )', [omega[i], osc_strength[i]])

plt.tight_layout()
plt.savefig('spectrumHHG.png', transparent=transparency)
st.pyplot(fig)
with open('spectrumHHG.png', 'rb') as f:
     st.download_button('Download plot as PNG file', f, file_name='spectrumHHG.png')  

## Create a dataframe object with the spectrum data
chart_data = pd.DataFrame(omega, columns=['Energy ('+energy_units+')'])
chart_data.insert(loc=1, column='Intensity', value=osc_strength)
# chart_data = chart_data.set_index('Frequency')
# st.line_chart(chart_data)
st.write('### OUTPUT RT-TDDFT HHG Spectrum (hhgspec)')
st.write(chart_data)
# determining the name of the file
export_rtspec_file_name = 'hhgspec.xlsx'
# saving to excel file
chart_data.to_excel(export_rtspec_file_name)
with open(export_rtspec_file_name, 'rb') as f:
     st.download_button('Download EXCEL file', f, file_name=export_rtspec_file_name)  
# st.write('rtspec is written to Excel File successfully.')
rtspec_str = chart_data.to_string(header=False, index=False)
rtspec_str = '$hhgspec\nEnergy ('+energy_units+')\t Intensity\n' + rtspec_str
rtspec_str = rtspec_str + '\n$end'

st.write('### OUTPUT (hhgspec)')
st.text_area(label='Absorption Spectra in `rtspec` format', value = rtspec_str, height=400)
