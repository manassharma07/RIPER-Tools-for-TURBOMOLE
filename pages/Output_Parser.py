import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO

def parse_energies(text):
    kinetic_energy = []
    coulomb_energy = []
    exchange_corr_energy = []
    total_energy = []
    
    lines = text.split('\n')
    for line in lines:
        if "KINETIC ENERGY" in line:
            kinetic_energy.append(float(line.split()[4]))
        elif "COULOMB ENERGY" in line:
            coulomb_energy.append(float(line.split()[4]))
        elif "EXCH. & CORR. ENERGY" in line:
            exchange_corr_energy.append(float(line.split()[4]))
        elif "TOTAL ENERGY" in line:
            total_energy.append(float(line.split()[4]))
    
    return kinetic_energy, coulomb_energy, exchange_corr_energy, total_energy


st.title("RIPER Output Parser")

st.write('You can either paste the output file contents below or upload the source file')
contents = st.text_area(label='Enter the contents of the output file here', value='', placeholder='Put your text here',
                        height=400, key='input_text_area')
# Create a file uploader widget
file = st.file_uploader("or Upload the file")

if file is not None:
    # If a file is uploaded, read its contents
    # contents = file.read()
    # To read file as bytes:
    bytes_data = file.getvalue()

    # To convert to a string based IO:
    stringio = StringIO(file.getvalue().decode("utf-8"))

    # To read file as string:
    contents = stringio.read()

if contents != '':
    file_contents = contents #upload_file.read()
    energies = parse_energies(file_contents)
    
    st.subheader("Parsed Energies")
    data = {
        "SCF Iteration": list(range(1, len(energies[0]) + 1)),
        "Kinetic Energy": energies[0],
        "Coulomb Energy": energies[1],
        "Exchange Corr. Energy": energies[2],
        "Total Energy": energies[3]
    }
    df = pd.DataFrame(data)
    st.dataframe(df)

    st.subheader("Energy vs SCF Iteration")
    plt.figure(figsize=(10, 6))
    plt.plot(data["SCF Iteration"], data["Total Energy"], marker='o', linestyle='-', color='b', label='Total Energy')
    plt.xlabel("SCF Iteration")
    plt.ylabel("Energy")
    plt.title("Energy vs SCF Iteration")
    plt.legend()
    st.pyplot()


