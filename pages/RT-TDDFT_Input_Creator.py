import streamlit as st

def generate_field_input(selected_field, amplitude_x, amplitude_y, amplitude_z, tzero, width, omega, sigma, phase_x, phase_y, phase_z):
    if selected_field == "Static":
        return f"$fields\n    electric on\n$electric field\n    amplitude x={amplitude_x}  y={amplitude_y}  z={amplitude_z}\n    static"
    elif selected_field == "Gaussian":
        return f"$fields\n    electric on\n$electric field\n    amplitude x={amplitude_x}  y={amplitude_y}  z={amplitude_z}\n    gaussian  tzero={tzero}  width={width}"
    elif selected_field == "Laser":
        return f"$fields\n    electric on\n$electric field\n    amplitude x={amplitude_x}  y={amplitude_y}  z={amplitude_z}\n    phase x={phase_x}  y={phase_y}  z={phase_z}\n  laser  omega={omega}  sigma={sigma}"

def generate_rttddft_input(magnus, scf, iterlim, time, tstep, print_step, damping, min_energy, max_energy, energy_step, print_density, print_energy, print_dipole, selected_field):
    rttddft_input = f"$rttddft\n    magnus {magnus}\n    scf {scf}\n    time {time}d0\n    tstep {tstep}d0\n    print step {print_step}\n"
    if scf=="on":
        rttddft_input += f"    iterlim {iterlim}"
    rttddft_input += f"    damping {damping}d0\n    min energy = {min_energy}d0\n    max energy = {max_energy}d0\n    energy step {energy_step}d0"
    if print_density:
        rttddft_input += f"$rtdens\n$pointvalper  fmt=cub\n    dens"
    if print_energy:
        rttddft_input += f"$rtenergy"
    if print_dipole:
        rttddft_input += f"$rtdipol"
    if selected_field=="Gaussian":
        rttddft_input += f"$rtspec"
    return rttddft_input


st.title("Input Creation for RT-TDDFT Simulation")

st.write('## Set Electric Field Parameters')

field_types = ["Static", "Gaussian", "Laser"]
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
    st.checkbox('Calculate and save absorption sepctrum to `rtspec` file? (only possible for Gaussian field)')

col1_rt, col2_rt, col3_rt = st.columns(3)
print_energy = col1_rt.checkbox('Print energy (`rtenrgy`) at each time step for post-processing?', value=True)
print_dipole = col2_rt.checkbox('Print dipole moment (`rtdipo`) at each time step for post-processing?', value=True)
print_density = col3_rt.checkbox('Print density at each time step for post-processing? (Will take up some disk space)', value=True)

magnus = col1_rt.selectbox("Magnus Expansion Order", [2, 4], index=0)
scf = col2_rt.radio("Use SCF Procedure for Time-Integration?", ["on", "off #(Use Predictor Corrector Scheme)"], index=1)
if scf=="on":
    iterlim = col3_rt.number_input("Max SCF Cycles", value=15)
time = col1_rt.number_input("Evolution Time (au)", value=1000.0)
tstep = col2_rt.number_input("Time Step (au)", value=0.1)
print_step = col3_rt.number_input("Print Step", value=1)
min_energy = col1_rt.number_input("Min Energy (au)", value=0.00)
max_energy = col2_rt.number_input("Max Energy (au)", value=0.75)
energy_step = col3_rt.number_input("Energy Step (au)", value=0.005)
damping = col1_rt.number_input("Damping Factor", value=0.004)

st.write("## Generated RT-TDDFT Input:")
rttddft_input = generate_rttddft_input(magnus, scf, iterlim, time, tstep, print_step, damping, min_energy, max_energy, energy_step, print_density, print_energy, print_dipole, selected_field, selected_field)
st.text_area(label="Add the following to the `control` file:", value=rttddft_input, height=400)