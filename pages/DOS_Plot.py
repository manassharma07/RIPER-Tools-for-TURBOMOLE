import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt


# Set page config
st.set_page_config(page_title='Convergence Tips for RIPER', layout='wide', page_icon="⚛️",
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


# Title
st.title("Density of States Plotter")

# Snippet of code for TURBOMOLE
st.subheader("RIPER `control` file snippet:")
st.code("""
# Add this snippet to your control file
$dosper width=real emin=real emax=real scal=real npt=integer
""")
st.write('and run `riper` as')
st.code('nohup riper -proper > dos.out &', language='shell')

# Input for spin-restricted or spin-unrestricted calculation
spin_restricted = st.checkbox("Spin-Restricted Calculation", value=True)

# Input for DOS file contents
st.subheader("Paste DOS File Contents:")

if spin_restricted:
    dos_contents = st.text_area("DOS File Contents", height=400)
else:
    dos_col1, dos_col2 = st.columns(2)
    dos_contents_alpha = dos_col1.text_area("DOS File (dos_alpha) Contents", height=400)
    dos_contents_beta = dos_col2.text_area("DOS File (dos_beta) Contents", height=400)


# Convert DOS contents to DataFrame
if (spin_restricted and dos_contents) or (not spin_restricted and dos_contents_alpha and dos_contents_beta):

    if spin_restricted:
        dos_lines = dos_contents.strip().split('\n')
        header_line = dos_lines[1].strip().split()[1:]
        data_lines = [line.strip().split() for line in dos_lines[2:]]
    else:
        dos_lines_alpha = dos_contents_alpha.strip().split('\n')
        header_line_alpha = dos_lines_alpha[1].strip().split()[1:]
        data_lines_alpha = [line.strip().split() for line in dos_lines_alpha[2:]]

        dos_lines_beta = dos_contents_beta.strip().split('\n')
        header_line_beta = dos_lines_beta[1].strip().split()[1:]
        data_lines_beta = [line.strip().split() for line in dos_lines_beta[2:]]

    if spin_restricted:
        dos_df = pd.DataFrame(data_lines, columns=header_line)
        dos_df = dos_df.astype(float)
    else:
        dos_df_alpha = pd.DataFrame(data_lines_alpha, columns=header_line_alpha)
        dos_df_alpha = dos_df_alpha.astype(float)

        dos_df_beta = pd.DataFrame(data_lines_beta, columns=header_line_beta)
        dos_df_beta = dos_df_beta.astype(float)

    if spin_restricted:
        # Convert energy to eV and shift by Fermi level
        fermi_level_au = st.number_input("Fermi Level (atomic units)", step=0.01)
        fermi_level_ev = fermi_level_au * 27.211324570273
        dos_df['energy'] = (dos_df['energy'] - fermi_level_au) * 27.211324570273
    else:
        # Convert energy to eV and shift by Fermi level
        fermi_level_au_alpha = dos_col1.number_input("Fermi Level Alpha (atomic units)", step=0.01)
        fermi_level_ev_alpha = fermi_level_au_alpha * 27.211324570273
        dos_df_alpha['energy'] = (dos_df_alpha['energy'] - fermi_level_au_alpha) * 27.211324570273
        # Convert energy to eV and shift by Fermi level
        fermi_level_au_beta = dos_col2.number_input("Fermi Level Beta (atomic units)", step=0.01)
        fermi_level_ev_beta = fermi_level_au_beta * 27.211324570273
        dos_df_beta['energy'] = (dos_df_beta['energy'] - fermi_level_au_beta) * 27.211324570273


    st.write('### Parsed DOS with energy in eV and shifted so that the Fermi energy is at 0 eV')
    if spin_restricted:
        st.dataframe(dos_df, height=400)
    else:
        dos_col1.dataframe(dos_df_alpha, height=400)
        dos_col2.dataframe(dos_df_beta, height=400)

    if spin_restricted:
        # saving to excel file
        excel_filename = 'dos.xlsx'
        dos_df.to_excel(excel_filename)
        with open(excel_filename, 'rb') as f:
            st.download_button('Download EXCEL file', f, file_name=excel_filename)  
    else:
        # saving to excel file
        excel_filename = 'dos_alpha.xlsx'
        dos_df_alpha.to_excel(excel_filename)
        with open(excel_filename, 'rb') as f:
            dos_col1.download_button('Download EXCEL file', f, file_name=excel_filename) 
        # saving to excel file
        excel_filename = 'dos_beta.xlsx'
        dos_df_beta.to_excel(excel_filename)
        with open(excel_filename, 'rb') as f:
            dos_col2.download_button('Download EXCEL file', f, file_name=excel_filename) 

    # Plotting
    st.subheader("DOS Plot:")

    col1, col2 = st.columns(2)
    transparency = col1.checkbox('Make the plot transparent', value=True)
    line_width = col2.slider('Line Width', 0.4, 10., 2.1, step=0.1)
    isGrids = col1.checkbox('Grids', value=False)
    title = col1.text_input("Title", "Density of States")
    fontsize_title = col2.number_input("Title Font Size", value=18)
    st.divider()

    col1, col2 = st.columns(2)
    if spin_restricted:
        xmin = col1.number_input("X-axis Min", value=dos_df['energy'].min())
        xmax = col2.number_input("X-axis Max", value=dos_df['energy'].max())

        ymin = col1.number_input("Y-axis Min", value=dos_df['total'].min())
        ymax = col2.number_input("Y-axis Max", value=dos_df['total'].max())
    else:
        xmin = col1.number_input("X-axis Min", value=min(dos_df_alpha['energy'].min(), dos_df_beta['energy'].min()))
        xmax = col2.number_input("X-axis Max", value=max(dos_df_alpha['energy'].max(), dos_df_beta['energy'].max()))

        ymin = col1.number_input("Y-axis Min", value=-dos_df_beta['total'].max())
        ymax = col2.number_input("Y-axis Max", value=dos_df_alpha['total'].max())

    xaxis_title = col1.text_input("X-axis Title", "Energy (eV)")
    yaxis_title = col2.text_input("Y-axis Title", "DOS (states/eV/unit cell)")

    fontsize_xlabel = col1.number_input("X-axis Title Font Size", value=14)
    fontsize_ylabel = col2.number_input("Y-axis Title Font Size", value=14)

    fontsize_xticks = col1.number_input("X-axis Tick Labels Font Size", value=14)
    fontsize_yticks = col2.number_input("Y-axis Tick Labels Font Size", value=12)

    figsize_x = col1.number_input("Plot Size (X)", value=10)
    if spin_restricted:
        figsize_y = col2.number_input("Plot Size (Y)", value=4)
    else:
        figsize_y = col2.number_input("Plot Size (Y)", value=5)

    fontsize_legend = col1.number_input("Legend Font Size", value=14)
    draw_fermi = col2.checkbox('Draw Vertical Line at Fermi Level', value=True)



    plt.figure(figsize=(figsize_x, figsize_y))
    if isGrids:
        plt.grid(True)
    else:
        plt.grid(False)

    if spin_restricted:
        plt.fill_between(dos_df['energy'], dos_df['total'], color='lightgrey', label='Total DOS')

        num_columns = len(header_line) - 1
        # colors = plt.cm.viridis.colors[:num_columns]
        colors = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494']
        for i, col in enumerate(header_line[2:]):
            plt.plot(dos_df['energy'], dos_df[col], color=colors[i], label=col, lw=line_width)
    else:
        plt.fill_between(dos_df_alpha['energy'], dos_df_alpha['total'], color='lightgrey', label='Total DOS')

        num_columns = len(header_line_alpha) - 1
        # colors = plt.cm.viridis.colors[:num_columns]
        colors = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494']
        for i, col in enumerate(header_line_alpha[2:]):
            plt.plot(dos_df_alpha['energy'], dos_df_alpha[col], color=colors[i], label=col, lw=line_width)

        plt.fill_between(dos_df_beta['energy'], -dos_df_beta['total'], color='lightgrey')

        num_columns = len(header_line_beta) - 1
        # colors = plt.cm.viridis.colors[:num_columns]
        colors = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494']
        for i, col in enumerate(header_line_beta[2:]):
            plt.plot(dos_df_beta['energy'], -dos_df_beta[col], color=colors[i], lw=line_width)

    
    plt.axvline(x=0, color='black', linestyle='--', label='Fermi Level', lw=0.4)

    plt.xlabel(xaxis_title, fontsize=fontsize_xlabel)
    plt.ylabel(yaxis_title, fontsize=fontsize_ylabel)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    

    plt.legend(fontsize=fontsize_legend)
    plt.title(title, fontsize=fontsize_title)

    # Customize labels and ticks

    plt.xticks(fontsize=fontsize_xticks)
    plt.yticks(fontsize=fontsize_yticks)

    plt.savefig("dos_plot.png", format='png', transparent=transparency, dpi=400)
    st.pyplot(plt)

    # Download plot as PNG
    st.subheader("Download Plot:")
    with open('dos_plot.png', 'rb') as f:
        st.download_button('Download plot as PNG file', f, file_name='dos_plot.png')
