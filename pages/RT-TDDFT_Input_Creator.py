import streamlit as st

def generate_input(selected_field, amplitude_x, amplitude_y, amplitude_z, tzero, width, omega, sigma, phase_x, phase_y, phase_z):
    if selected_field == "Static":
        return f"$fields\nelectric on\n$electric field\namplitude x={amplitude_x} y={amplitude_y} z={amplitude_z}\nstatic"
    elif selected_field == "Gaussian":
        return f"$fields\nelectric on\n$electric field\namplitude x={amplitude_x} y={amplitude_y} z={amplitude_z}\ngaussian tzero={tzero} width={width}"
    elif selected_field == "Laser":
        return f"$fields\nelectric on\n$electric field\namplitude x={amplitude_x} y={amplitude_y} z={amplitude_z}\nphase x={phase_x} y={phase_y} z={phase_z}\nlaser omega={omega} sigma={sigma}"




st.title("Input Creation for RT-TDDFT Simulation")

field_types = ["Static", "Gaussian", "Laser"]
selected_field = st.selectbox("Select Field Type", field_types)

st.write("Enter Electric Field Parameters:")

col1, col2, col3 = st.columns(3)
amplitude_x = col1.text_input(label = 'Amplitude along x (a.u.)', value='2.0E-5',key='Ex')
amplitude_y = col2.text_input(label = 'Amplitude along y (a.u.)', value='2.0E-5',key='Ey')
amplitude_z = col3.text_input(label = 'Amplitude along z (a.u.)', value='2.0E-5',key='Ez')

if selected_field == "Gaussian":
    col1_gaussian, col2_gaussian = st.columns(2)
    tzero = col1_gaussian.number_input("Peak Position (tzero)", value=0.0)
    width = col2_gaussian.number_input("Peak Width (width)", value=0.0)

if selected_field == "Laser":
    omega = st.number_input("Frequency (omega)", value=0.0)
    sigma = st.number_input("FWHM (sigma)", value=0.0)
    phase_x = col1.number_input("Phase x", value=0.0)
    phase_y = col2.number_input("Phase y", value=0.0)
    phase_z = col3.number_input("Phase z", value=0.0)

st.write("Generated Electric Field Input:")
st.code(generate_input(selected_field, amplitude_x, amplitude_y, amplitude_z, tzero, width, omega, sigma, phase_x, phase_y, phase_z))

