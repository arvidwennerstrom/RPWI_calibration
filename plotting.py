# HK-plotting script
# By: Arvid Wennerström"



# ================================================
#               Imports & Stuff 
# ================================================
import numpy as np
from rpwi_data import rpwi_data  
from rpwi_data import calibration_coefficients
from math import ceil

# ================================================
#               Classes & Functions 
# ================================================
class Plot_options:
    def __init__(self,unit = 'TM',saturation_filter = True,align_Efield = True):
        self.unit = unit
        self.saturation_filter = saturation_filter
        self.align_Efield = align_Efield
        

def tt2000_to_readable(tt2000,precision = 9):
    # Convert TT2000 (nanoseconds) to ephemeris time (TDB) in seconds
    ephem_time = spice.unitim(tt2000/ 1e9, 'TT', 'TDB')
    
    # Convert ephemeris time to human-readable UTC string
    utcs = spice.spiceypy.et2utc(ephem_time , format_str='ISOC', prec=precision)

    # # Extract only the time part (remove the date)
    # readable = utcs.split('T')[1]

    # Return all but year
    readable = utcs[5:]
    
    return readable


# Custom formatter function to convert displayed tt2000 values dynamically
def dynamic_time_formatter(x, pos):
    # Get the current limits of the x-axis to adjust precision based on zoom level
    current_xlim = plt.gca().get_xlim()
    range_x = current_xlim[1] - current_xlim[0]

    # Adjust precision based on zoom level
    if range_x > 1e12:  # Very zoomed out, show only date
        precision = 0
    elif range_x > 1e9:  # Show hours
        precision = 3
    elif range_x > 1e6:  # Show hours and minutes
        precision = 6
    else:  # Show full precision down to seconds
        precision = 9


    # Convert the tt2000 time to human-readable format with the determined precision
    readable_time = tt2000_to_readable(x,precision)
    return readable_time



# ================================================
#                   Main
# ================================================


# Input overall directory path, where folders 
# "datasets", "data_created" and "spice" exist
# ================================================ 
rootDir = "C:/Users/arvwe/Onedrive - KTH/MEX/IRF" # My laptop
# rootDir = "C:/Users/arvidwen/Onedrive - KTH/MEX/IRF" # KTH computers
date = '240823'


#           Choose settings for plot
# ================================================
plot_matplotlib = True 
plot_plotly = False

# unit: 'TM', 'U' (voltage), 'EF' (Electric field)
# saturation filte: True/False
# Align E-field: True/False
plot_options = Plot_options('TM',False,False)

mask_threshold = 3



#               Load data
# ================================================
loaded_LP = np.load(rootDir + "/data_created/LP-SID1_" + date + ".npz")
Epoch = loaded_LP["Epoch"]
LP_data = loaded_LP["LP_data"]
Mask = loaded_LP["Mask"]



#   Add SPICE-kernels, for converting time to human-readable
# ================================================
import spiceypy as spice    
ls_spice_kernels = [
    (rootDir + '/spice/JUICE/kernels/lsk/naif0012.tls'),
    (rootDir + '/spice/JUICE/kernels/sclk/juice_step_240828_v01.tsc')]
spice.furnsh(ls_spice_kernels)



#     Transform E-field to SC x-,y- and z-axes
# ================================================
if plot_options.align_Efield:
    # M is transformation matrix in: U = M*E, where U is voltage and E is electric field
    
    # NOTE: This is for mode 0
    M = [rpwi_data.LP12_distance,
         rpwi_data.LP23_distance,
         rpwi_data.LP34_distance]
    M_inv = np.linalg.inv(M)

    E_field = np.matmul(M_inv,LP_data[:3])


#          Scale and update y-axis values
# ================================================
if plot_options.unit == 'U' or plot_options.unit == 'EF':
    LP_data = calibration_coefficients*LP_data
    ylabel_diff = ['V', 'V']
    if plot_options.unit == 'EF':
        LP_data[0] = 1e3*LP_data[0]/np.linalg.norm(rpwi_data.LP12_distance)
        LP_data[1] = 1e3*LP_data[1]/np.linalg.norm(rpwi_data.LP23_distance)
        LP_data[2] = 1e3*LP_data[2]/np.linalg.norm(rpwi_data.LP34_distance)
        ylabel = ['mV/m', 'V']
else:
    ylabel = ['TM', 'TM']


if plot_matplotlib:
    #               Plot (using matplotlib)
    # ================================================
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    
    # Divide the data into chunks, by number of data points in "limit"
    limit = 5e6
    number_of_data_chunks = ceil(len(Epoch)/limit)
    chunk_size = round(len(Epoch)/number_of_data_chunks)




    #           Repeat for each chunk
    # ================================================
    for chunk_i in range(2):
        # Decide which data indecies shall be plotted in figure
        lower_i = int(chunk_size*chunk_i)
        upper_i = int(chunk_size*(chunk_i+1)-1)
            

        # Get the chunk of data to plot in this figure 
        # ================================================
        Epoch_chunk = Epoch[lower_i:upper_i]
        LP_data_chunk = LP_data[:,lower_i:upper_i]
        Mask_chunk = Mask[:,lower_i:upper_i]


        # Plot probe differentials and single ended probe 
        # ================================================
        fig, axes = plt.subplots(2,1, sharex=True)
        titles = ['Probe differentials','Single ended probe (P4)']
        labels = ['P12', 'P23', 'P34', 'P4']

        
        [axes[0].plot(Epoch_chunk[Mask_chunk[i] >= mask_threshold],LP_data_chunk[i][Mask_chunk[i] >= mask_threshold],'o', markersize = 1, linestyle = 'None', label=labels[i]) for i in range(3)]
        axes[1].plot(Epoch_chunk[Mask_chunk[3] >= mask_threshold],LP_data_chunk[3][Mask_chunk[3] >= mask_threshold],'o', markersize = 1, linestyle = 'None', label=labels[3])


        # Setup for each of the subplots
        for subplot_i in range(2):
            # Set custom formatter for x-axis
            axes[subplot_i].xaxis.set_major_formatter(FuncFormatter(dynamic_time_formatter))

            # Add labels and legend
            axes[subplot_i].set_title(titles[subplot_i])
            axes[subplot_i].set_xlabel('Time')
            axes[subplot_i].set_ylabel(ylabel[subplot_i])
            axes[subplot_i].legend(loc='upper right')

            # Show grid
            axes[subplot_i].grid(True)

        # Change window title
        window_title = '20' + date[0:2] + '-' + date[2:4] + '-' + date[4:6] + ' Chunk number ' + str(chunk_i+1)
        fig.canvas.manager.set_window_title(window_title)

        # Rotate the x-axis labels for better readability
        plt.gcf().autofmt_xdate()
        

        if plot_options.align_Efield:
            plt.figure()
            plt.plot(Epoch[lower_i:upper_i],E_field[0][lower_i:upper_i],'o', markersize = 1, linestyle = 'None', label='X')
            plt.plot(Epoch[lower_i:upper_i],E_field[1][lower_i:upper_i],'o', markersize = 1, linestyle = 'None', label='Y')
            plt.plot(Epoch[lower_i:upper_i],E_field[2][lower_i:upper_i],'o', markersize = 1, linestyle = 'None', label='Z')
            plt.title('SC axes aligned E-field')
            plt.xlabel('Time')
            plt.ylabel(ylabel[0])
            plt.legend(loc='upper right')
            plt.grid(True)
            # plt.canvas.manager.set_window_title(window_title)
            # Rotate the x-axis labels for better readability
            plt.gcf().autofmt_xdate()
    

    # Show all the plots
    # ================================================
    plt.show()



if plot_plotly:
    #               Plot (using plotly)
    # ================================================
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots

    fig = make_subplots(rows=1, cols=1, shared_xaxes=True)

    # Create plotly express scatter plots
    scatter_1 = px.scatter(x=Epoch, y=LP_data[0], labels={'y': 'P12'})
    scatter_2 = px.scatter(x=Epoch, y=LP_data[1], labels={'y': 'P23'})
    # scatter_3 = px.scatter(x=Epoch[lower_i:upper_i], y=LP_data[2][lower_i:upper_i], labels={'y': 'P34'})
    # scatter_4 = px.scatter(x=Epoch[lower_i:upper_i], y=LP_data[3][lower_i:upper_i], labels={'y': 'P4'})

    # Add scatter plots to subplots
    fig.add_trace(scatter_1.data[0], row=1, col=1)
    fig.add_trace(scatter_2.data[0], row=1, col=1)
    # fig.add_trace(scatter_3.data[0], row=1, col=1)
    # fig.add_trace(scatter_4.data[0], row=2, col=1)

    # Update layout for the figure
    fig.update_layout(title_text="LP Data with Shared X-Axis")

    # Show grid for both subplots
    fig.update_xaxes(showgrid=True)
    fig.update_yaxes(showgrid=True)

    # Add x-axis and y-axis labels
    fig.update_xaxes(title_text="Time", row=2, col=1)
    fig.update_yaxes(title_text=ylabel[0], row=1, col=1)
    fig.update_yaxes(title_text=ylabel[1], row=2, col=1)

    # Show the plot
    fig.show()


 
    # # Melt the DataFrame to a long format
    # df_melted = df.melt(id_vars='x', value_vars=['y1', 'y2', 'y3', 'y4'], var_name='Series', value_name='y')


    # # Create scatter plot
    # fig = px.scatter(df_melted, x='x', y ='y', color='Series')

    # # Show the plot
    # fig.show()

