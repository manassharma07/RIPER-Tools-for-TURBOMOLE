import streamlit as st

def generate_input(selected_field, amplitude_x, amplitude_y, amplitude_z, tzero, width, omega, sigma, phase_x, phase_y, phase_z):
    if selected_field == "Static":
        return f"$fields\n  electric on\n$electric field\n  amplitude x={amplitude_x}  y={amplitude_y}  z={amplitude_z}\n  static"
    elif selected_field == "Gaussian":
        return f"$fields\n  electric on\n$electric field\n  amplitude x={amplitude_x}  y={amplitude_y}  z={amplitude_z}\n  gaussian  tzero={tzero}  width={width}"
    elif selected_field == "Laser":
        return f"$fields\n  electric on\n$electric field\n  amplitude x={amplitude_x}  y={amplitude_y}  z={amplitude_z}\n  phase x={phase_x}  y={phase_y}  z={phase_z}\n  laser  omega={omega}  sigma={sigma}"




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

st.write('## Input Text')
st.text_area(label="Add the following to the `control` file:", value=generate_input(selected_field, amplitude_x, amplitude_y, amplitude_z, tzero, width, omega, sigma, phase_x, phase_y, phase_z), height=200)

