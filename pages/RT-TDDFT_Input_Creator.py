import streamlit as st
import pandas as pd

def generate_field_input(selected_field, amplitude_x, amplitude_y, amplitude_z, tzero, width, omega, sigma, phase_x, phase_y, phase_z):
    if selected_field == "Static":
        return f"$fields\n    electric on\n$electric field\n    amplitude x={amplitude_x}  y={amplitude_y}  z={amplitude_z}\n    static"
    elif selected_field == "Gaussian":
        return f"$fields\n    electric on\n$electric field\n    amplitude x={amplitude_x}  y={amplitude_y}  z={amplitude_z}\n    gaussian  tzero={tzero}  width={width}"
    elif selected_field == "Laser":
        return f"$fields\n    electric on\n$electric field\n    amplitude x={amplitude_x}  y={amplitude_y}  z={amplitude_z}\n    phase x={phase_x}  y={phase_y}  z={phase_z}\n  laser  omega={omega}  sigma={sigma}"

def generate_rttddft_input(magnus, scf, iterlim, time, tstep, print_step, damping, min_energy, max_energy, energy_step, print_density, print_energy, print_dipole, selected_field, print_spectrum):
    rttddft_input = f"$rttddft\n    magnus {magnus}\n    scf {scf}\n    time {time}d0\n    tstep {tstep}d0\n    print step {print_step}\n"
    if scf=="on":
        rttddft_input += f"    iterlim {iterlim}"
    rttddft_input += f"    damping {damping}d0\n    min energy = {min_energy}d0\n    max energy = {max_energy}d0\n    energy step {energy_step}d0"
    if print_density:
        rttddft_input += f"\n$rtdens\n$pointvalper  fmt=cub\n    dens"
    if print_energy:
        rttddft_input += f"\n$rtenergy"
    if print_dipole:
        rttddft_input += f"\n$rtdipol"
    if selected_field=="Gaussian" and print_spectrum:
        rttddft_input += f"\n$rtspectrum  eV"
    return rttddft_input


st.title("Input Creation for RT-TDDFT Simulation")

st.write('Using this tool you can create input for running an RT-TDDFT calculation using the `RIPER` module.')
st.write('First you need to specify some parameters for the Electric Field and paste the resulting output to your `control` file.')
st.write('Next, you need to specify the RT-TDDFT parameters and paste the resulting output to the `control` file.')
st.warning('Please Note: The current defaults when you load the page are suitable for performing an RT-TDDFT calculation with the goal of calculating the absorption spectrum. The Gaussian field perturbs the molecule in all three directions and the density is evolved for 24.18 fs.')

st.write('## Set Electric Field Parameters')

field_types = ["Gaussian", "Laser", "Static"]
selected_field = st.selectbox("Select Field Type", field_types)

st.write("Enter Electric Field Parameters:")
tzero = None
width = None
sigma = None
omega = None
phase_x = None
phase_y = None
phase_z = None

col1, col2, col3 = st.columns(3)
amplitude_x = col1.text_input(label = 'Amplitude along x (a.u.)', value='2.0E-5',key='Ex')
amplitude_y = col2.text_input(label = 'Amplitude along y (a.u.)', value='2.0E-5',key='Ey')
amplitude_z = col3.text_input(label = 'Amplitude along z (a.u.)', value='2.0E-5',key='Ez')

if selected_field == "Gaussian":
    col1_gaussian, col2_gaussian = st.columns(2)
    tzero = col1_gaussian.text_input("Peak Position (tzero in a.u.)", value='3.0')
    width = col2_gaussian.text_input("Peak Width (width in a.u.)", value='0.2')

if selected_field == "Laser":
    col1_laser, col2_laser = st.columns(2)
    omega = col1_laser.text_input("Frequency/a.u. (omega)", value='0.04556916') # 1.24 eV
    sigma = col2_laser.text_input("FWHM/a.u. (sigma)", value='1379.0')
    phase_x = col1.text_input("Phase x (radians)", value='0.0')
    phase_y = col2.text_input("Phase y (radians)", value='0.0')
    phase_z = col3.text_input("Phase z (radians)", value='0.0')

st.write('## Input Text for Electric Field')
st.text_area(label="Add the following to the `control` file:", value=generate_field_input(selected_field, amplitude_x, amplitude_y, amplitude_z, tzero, width, omega, sigma, phase_x, phase_y, phase_z), height=200)

st.write('## RT-TDDFT Options')

if selected_field=='Gaussian':
    print_spectrum = st.checkbox('Calculate and save absorption sepctrum to `rtspec` file? (only possible for Gaussian field)', value=True)
else:
    print_spectrum = False

col1_rt, col2_rt, col3_rt = st.columns(3)
print_energy = col1_rt.checkbox('Print energy (`rtenrgy`) at each time step for post-processing?', value=True)
print_dipole = col2_rt.checkbox('Print dipole moment (`rtdipo`) at each time step for post-processing?', value=True)
print_density = col3_rt.checkbox('Print density at each time step for post-processing? (Will take up some disk space)', value=False)

magnus = col1_rt.selectbox("Magnus Expansion Order", [2, 4], index=0)
scf = col2_rt.radio("Use SCF Procedure for Time-Integration?", ["on", "off #(Use Predictor Corrector Scheme)"], index=1)
if scf=="on":
    iterlim = col3_rt.number_input("Max SCF Cycles", value=15)
else:
    iterlim = None
time = col1_rt.number_input("Evolution Time (au)", value=1000.0)
tstep = col2_rt.number_input("Time Step (au)", value=0.1)
print_step = col3_rt.number_input("Print Step", value=1)
min_energy = col1_rt.number_input("Min Energy (au)", value=0.00)
max_energy = col2_rt.number_input("Max Energy (au)", value=0.75)
energy_step = col3_rt.number_input("Energy Step (au)", value=0.005)
damping = col1_rt.number_input("Damping Factor", value=0.004)

st.write("## Generated RT-TDDFT Input:")
rttddft_input = generate_rttddft_input(magnus, scf, iterlim, time, tstep, print_step, damping, min_energy, max_energy, energy_step, print_density, print_energy, print_dipole, selected_field, print_spectrum)
st.text_area(label="Add the following to the `control` file:", value=rttddft_input, height=400)


# Provide some hints
hints_data = [
    "The absorption spectrum can only be calculated with the Gaussian Field.",
    "The defaults when you launch this web app should be suitable for calculating the absorption spectrum of a molecule. The spectrum is saved to the `rtspec` file in the calculation directory.",
    "In some cases you may need to increase the amplitudes of the Gaussian pulse very slightly.",
    "A higher evolution time means a sharper absorption spectrum with more pronounced peaks.",
    "You should always use the `$rtenergy` keyword as it helps you monitor your RT-TDDFT simulation more closely. You can use it to see if the energies are diverging, which can happen when using the Predictor-Corrector scheme for time-integration. ",
    "Predictor-Corrector (PC) scheme (`scf off`) requires less KS matrix builds and should be faster than SCF scheme (`scf on`). However, for a larger time-step PC can become unstable.",
    "In tests it was found that the time step for RT-TDDFT could be safely increased to 0.5 a.u. to accelerate simulations with the PC scheme. But you will need to test this on a case-by-case basis. The default 0.1 a.u. is a bit conservative.",
    "In tests it was found that a time step that is 20% of the theoretical maximum time step (\pi / \omega_{max}) should be a safe choice. Here, \omega_{max} is the maximum frequency upto which the absorption spectrum would be plotted."
]

# Create a DataFrame from the hints data
hints_df = pd.DataFrame(hints_data, columns=["Hints"])

st.write("## Useful Hints")

# Display the hints table
st.table(hints_df)

st.write('## Keywords and their Meanings')
keywords = {
    "Keyword": ["magnus", "scf", "iterlim", "time", "tstep", "print step", "damping", "min energy", "max energy", "energy step"],
    "Meaning": [
        "Can take values 2 or 4. '2' for second order Magnus expansion and '4' for fourth order Magnus expansion. Default value is 2 and good enough for most applications.",
        "If `on`, then SCF procedure is used for the time integration. If off then Predictor-Corrector scheme is used instead.",
        "Max SCF cycles if scf is `on`. Default value is 15.",
        "Specifies the evolution time in au. (1 au= 0.02419 fs)",
        "The time step for the time evolution in au. 0.1 au is usually a good starting point.",
        "Specifies the number of steps n after which the dipole moments and energies are printed out if requested. Default value is 100. That means the quantities are printed out at every 100 steps. To have all the information for post-processing, a value of 1 is recommended.",
        "Only valid for absorption spectrum calculation. It is the factor gamma in the equation to calculate the complex polarizability tensor. Default value is 0.004 au. Recommended values in the range of 0.003 au to 0.005 au.",
        "Only valid for absorption spectrum calculation. Specifies the minimum value of the energy range used to perform the Fourier transform from time to frequency space. Units: au. Default value is 0.15 au.",
        "Only valid for absorption spectrum calculation. Specifies the maximum value of the energy range used to perform the Fourier transform from time to frequency space. Units: au. Default value is 0.625 au.",
        "Only valid for absorption spectrum calculation. Specifies the step value or energy interval dE at which to sample the energy values for Fourier transform and absorption spectrum plotting. Units: au. Default value is 0.005 au."
    ]
}

keywords_df = pd.DataFrame(keywords)
# Display the keywords table
st.table(keywords_df)

st.warning('## Implementation Paper')
st.write('Please cite the following paper if you use the RT-TDDFT feature')
st.write(f"Carolin Müller, Manas Sharma, Marek Sierka,  \n*Real-time time-dependent density functional theory using density fitting and the continuous fast multipole method*,  \n J Comput Chem. 2020; 41: 2573–258")

# Create a hyperlink for the citation
citation_link = "[Read the paper](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.26412)"

# Display the citation with the hyperlink
st.markdown(citation_link, unsafe_allow_html=True)