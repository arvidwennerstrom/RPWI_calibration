# RPWI Electric field data
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
    def __init__(self,style='-',unit = 'U',probe_differentials = False,potential_at_probes = False,Efield = False,Efield_xyz = False):
        self.style = style
        self.unit = unit
        self.probe_differentials = probe_differentials
        self.potential_at_probes = potential_at_probes
        self.Efield = Efield
        self.Efield_xyz = Efield_xyz
        

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
rootDir = "C:/Users/arvidwen/Onedrive - KTH/MEX/IRF" # KTH computers
date = '240820'


#           Choose settings for plot
# ================================================
plot_matplotlib = True 
plot_plotly = False

# Plot style:                   '-' (lines) or 'o' (dots)
# Unit:                         'TM'/'V' (voltage)
# Probe differentials (TM or V):     True/False
# SC potential at each probe:   True/False
# Electric field:               True/False
# SC axes aligned E-field:      True/False
plot_options = Plot_options('-','V',False,False,True,True)


# 1: TM saturation (either high or low).
# 2: Interference with RIME instrument.
# 3: Lunar wake data, will overwrite '1'.
# 4: High quality data, this does not contain any contaminations.
mask_limit = 4



#               Load data
# ================================================
loaded_LP = np.load(rootDir + "/data_created/LP-SID1_" + date + ".npz")
Epoch = loaded_LP["Epoch"]
LP_diffs_TM = loaded_LP["LP_data"]
Mask = loaded_LP["Mask"]

# NOTE: LP_data is now always converted to voltage from TM, so unit selection 
# above does nothing. This might be desirable to change in the future.
LP_diffs = calibration_coefficients*LP_diffs_TM



#   Add SPICE-kernels, for converting time to human-readable
# ================================================
import spiceypy as spice    
ls_spice_kernels = [
    (rootDir + '/spice/JUICE/kernels/lsk/naif0012.tls'),
    (rootDir + '/spice/JUICE/kernels/sclk/juice_step_240828_v01.tsc')]
spice.furnsh(ls_spice_kernels)



# Get SC potential at each probe, by adding differentials.
# Use this to approximate DC offset coefficients and get
# calibrated potentials.
# ================================================
LP_potentials = np.zeros(LP_diffs.shape)
LP_potentials[3] = LP_diffs[3]
LP_potentials[2] = LP_potentials[3] + LP_diffs[2]
LP_potentials[1] = LP_potentials[2] + LP_diffs[1]
LP_potentials[0] = LP_potentials[1] + LP_diffs[0]

# ================================================
# NOTE: Simple way of doing it, but is problematic when LP[3]/LP[n] --> inf.
k1 = np.mean(LP_potentials[3]/LP_potentials[0])
k2 = np.mean(LP_potentials[3]/LP_potentials[1])
k3 = np.mean(LP_potentials[3]/LP_potentials[2])


# NOTE: Alternative way, this removes outliers from LP[3]/LP[n] --> inf,
# by filtering by standard deviation
if False:
    from scipy import stats
    std_dev_limit = 1

    k1_all = LP_potentials[3]/LP_potentials[0]
    k2_all = LP_potentials[3]/LP_potentials[1]
    k3_all = LP_potentials[3]/LP_potentials[2]

    z_score1 = np.abs(stats.zscore(k1_all))
    z_score2 = np.abs(stats.zscore(k2_all))
    z_score3 = np.abs(stats.zscore(k3_all))

    k1 = np.mean(k1_all[z_score1 < std_dev_limit])
    k2 = np.mean(k2_all[z_score2 < std_dev_limit])
    k3 = np.mean(k3_all[z_score3 < std_dev_limit])
# ================================================


LP_potentials_calibr = np.zeros(LP_potentials.shape)
LP_potentials_calibr[3] = LP_potentials[3]
LP_potentials_calibr[2] = k3*LP_potentials[2]
LP_potentials_calibr[1] = k2*LP_potentials[1]
LP_potentials_calibr[0] = k1*LP_potentials[0]



# Calculate new differentials using calibrated potentials
# ================================================
LP_diffs_calibr = np.zeros(LP_potentials_calibr.shape)
LP_diffs_calibr[3] = LP_potentials_calibr[3]
LP_diffs_calibr[2] = LP_potentials_calibr[2] - LP_potentials_calibr[3]
LP_diffs_calibr[1] = LP_potentials_calibr[1] - LP_potentials_calibr[2]
LP_diffs_calibr[0] = LP_potentials_calibr[0] - LP_potentials_calibr[1]



# Get E-field between probes
# ================================================
if plot_options.Efield: 
    # 1e3's are used to give EF in mV/m, rather than V/m

    # NOTE: This is for mode 0
    EF_in_LPs = np.zeros((3,len(Epoch)))
    EF_in_LPs[0] = 1e3*LP_diffs[0]/np.linalg.norm(rpwi_data.LP12_distance)
    EF_in_LPs[1] = 1e3*LP_diffs[1]/np.linalg.norm(rpwi_data.LP23_distance)
    EF_in_LPs[2] = 1e3*LP_diffs[2]/np.linalg.norm(rpwi_data.LP34_distance)
    
    EF_in_LPs_calibr = np.zeros((3,len(Epoch)))
    EF_in_LPs_calibr[0] = 1e3*LP_diffs_calibr[0]/np.linalg.norm(rpwi_data.LP12_distance)
    EF_in_LPs_calibr[1] = 1e3*LP_diffs_calibr[1]/np.linalg.norm(rpwi_data.LP23_distance)
    EF_in_LPs_calibr[2] = 1e3*LP_diffs_calibr[2]/np.linalg.norm(rpwi_data.LP34_distance)



# Calculate E-field in SC x-,y- and z-axes
# ================================================
if plot_options.Efield_xyz:

    # M is transformation matrix in: U = M*E, where U is voltage and E is electric field
    M = [rpwi_data.LP12_distance,
         rpwi_data.LP23_distance,
         rpwi_data.LP34_distance]
    M_inv = np.linalg.inv(M)
    EF_in_xyz = 1e3*np.matmul(M_inv,LP_diffs[:3])
    EF_in_xyz_calibr = 1e3*np.matmul(M_inv,LP_diffs_calibr[:3])




#               Plot (using matplotlib)
# ================================================
if plot_matplotlib:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    
    
    # Divide the data into chunks, by number of data points in "limit".
    # Then plot a new figure for each chunk.
    # ================================================================
    limit = 4e8
    number_of_data_chunks = ceil(len(Epoch)/limit)
    chunk_size = round(len(Epoch)/number_of_data_chunks)


    for chunk_i in range(1):
        # Decide and get data that belongs to this chunk.
        # ================================================
        lower_i = int(chunk_size*chunk_i)
        upper_i = int(chunk_size*(chunk_i+1)-1)
                
        Epoch_chunk = Epoch[lower_i:upper_i]
        LP_diffs_chunk = LP_diffs[:,lower_i:upper_i]
        Mask_chunk = Mask[:,lower_i:upper_i]
        
        
        
        # General plot setup
        # ================================================
        window_title = '20' + date[0:2] + '-' + date[2:4] + '-' + date[4:6] + ' Chunk ' + str(chunk_i+1)
        linewidth = 0.3
        xlabel = "Timestamp (UTC)"


        # Plot SC potential at each probe
        # ================================================
        if plot_options.potential_at_probes:
            SC_potentials_chunk = LP_potentials[:,lower_i:upper_i]
            LP_potentials_calibr_chunk = LP_potentials_calibr[:,lower_i:upper_i]

            legends = [
                ['P01 (P12+P23+P34+P04)', 'P02 (P23+P34+P04)', 'P03 (P34+P04)', 'P04'],
                ['k1*U1', 'k2*U2', "k3*U3", 'U4']
                ]
            
            titles = [
                'SC potential at each probe',
                'Calibrated potentials'
            ]


            fig, axes = plt.subplots(2,1,sharex=True)
            [axes[0].plot(Epoch_chunk[Mask_chunk[i] >= mask_limit],SC_potentials_chunk[i][Mask_chunk[i] >= mask_limit], lw=linewidth, label=legends[0][i]) for i in range(4)]
            [axes[1].plot(Epoch_chunk[Mask_chunk[i] >= mask_limit],LP_potentials_calibr_chunk[i][Mask_chunk[i] >= mask_limit], lw=linewidth, label=legends[1][i]) for i in range(4)]
            

            # Setup for the subplots
            for i in range(len(axes)):
                axes[i].xaxis.set_major_formatter(FuncFormatter(dynamic_time_formatter))
                axes[i].set_xlabel(xlabel)
                axes[i].set_ylabel(plot_options.unit)
                axes[i].set_title(titles[i])
                axes[i].legend(loc='upper right')
                axes[i].grid(True)  
            fig.canvas.manager.set_window_title("LP potentials " + window_title)


            # Rotate the x-axis labels for better readability
            plt.gcf().autofmt_xdate()
    


        # Plot probe differentials (either in voltages or TM)
        # ================================================
        if plot_options.probe_differentials:
            LP_diffs_calibr_chunk = LP_diffs_calibr[:,lower_i:upper_i]

            legends = [
                'P12', 'P23', 'P34', 'P4']
            titles = [
                'Probe differentials',
                'Calibrated probe diff.',
                'Single ended (P4)'
            ]


            fig, axes = plt.subplots(3,1,sharex=True)
            if plot_options.style == '-':
                [axes[0].plot(Epoch_chunk[Mask_chunk[i] >= mask_limit],LP_diffs_chunk[i][Mask_chunk[i] >= mask_limit], lw=linewidth, label=legends[i]) for i in range(3)]
                [axes[1].plot(Epoch_chunk[Mask_chunk[i] >= mask_limit],LP_diffs_calibr_chunk[i][Mask_chunk[i] >= mask_limit], lw=linewidth, label=legends[i]) for i in range(3)]          
                axes[2].plot(Epoch_chunk[Mask_chunk[3] >= mask_limit],LP_diffs_calibr_chunk[3][Mask_chunk[3] >= mask_limit], lw = linewidth, label=legends[3])
            # elif plot_options.style == 'o':
            #     [plt.plot(Epoch_chunk[Mask_chunk[i] >= mask_limit],LP_diffs_chunk[i][Mask_chunk[i] >= mask_limit],'o', markersize = 1, linestyle = 'None',label=legends[i]) for i in range(4)]
            
            # Setup for the subplots
            for i in range(len(axes)):
                axes[i].xaxis.set_major_formatter(FuncFormatter(dynamic_time_formatter))
                axes[i].set_xlabel(xlabel)
                axes[i].set_ylabel(plot_options.unit)
                axes[i].set_title(titles[i])
                axes[i].legend(loc='upper right')
                axes[i].grid(True)  
            fig.canvas.manager.set_window_title('LP differentials ' + window_title)

            # Rotate the x-axis labels for better readability
            plt.gcf().autofmt_xdate()



        # Plot electric field strenth in probe differentials
        # ================================================
        if plot_options.Efield:
            EF_in_LPs_chunk = EF_in_LPs[:,lower_i:upper_i]
            EF_in_LPs_calibr_chunk = EF_in_LPs_calibr[:,lower_i:upper_i]

            fig, axes = plt.subplots(3,1, sharex=True)
            titles = ['Probe differential E-field','Calibrated diff. E-field','Single ended probe (P4)']
            legends = ['P12', 'P23', 'P34','P4']
            ylabels = ['mV/m','mV/m','V']

            if plot_options.style == '-':
                [axes[0].plot(Epoch_chunk[Mask_chunk[i] >= mask_limit],EF_in_LPs_chunk[i][Mask_chunk[i] >= mask_limit], lw = linewidth, label=legends[i]) for i in range(3)]
                [axes[1].plot(Epoch_chunk[Mask_chunk[i] >= mask_limit],EF_in_LPs_calibr_chunk[i][Mask_chunk[i] >= mask_limit], lw = linewidth, label=legends[i]) for i in range(3)]
                axes[2].plot(Epoch_chunk[Mask_chunk[3] >= mask_limit],LP_diffs_chunk[3][Mask_chunk[3] >= mask_limit], lw = linewidth, label=legends[3])

            # elif plot_options.style == 'o':
            #     [axes[0].plot(Epoch_chunk[Mask_chunk[i] >= mask_limit],LP_diffs_chunk[i][Mask_chunk[i] >= mask_limit],'o', markersize = 1, linestyle = 'None',legends=labels[i]) for i in range(3)]
            #     axes[2].plot(Epoch_chunk[Mask_chunk[3] >= mask_limit],LP_diffs_chunk[3][Mask_chunk[3] >= mask_limit],'o', markersize = 1, linestyle = 'None',legends=labels[3])


            [axes[i].xaxis.set_major_formatter(FuncFormatter(dynamic_time_formatter)) for i in range(len(axes))]
            [axes[i].set_ylabel(ylabels[i]) for i in range(len(axes))]
            [axes[i].legend(loc='upper right') for i in range(len(axes))]
            [axes[i].set_title(titles[i]) for i in range(len(axes))]
            [axes[i].grid(True) for i in range(len(axes))]
            
            fig.canvas.manager.set_window_title("EF " + window_title)

            # Rotate the x-axis labels for better readability
            plt.gcf().autofmt_xdate()



        # Plot electric field strenth in x,y,z-axes
        # ================================================
        if plot_options.Efield_xyz:
            EF_in_xyz_chunk = EF_in_xyz[:,lower_i:upper_i]
            EF_in_xyz_calibr_chunk = EF_in_xyz_calibr[:,lower_i:upper_i]


            fig, axes = plt.subplots(2,1, sharex=True)
            titles = ['Calibrated axis aligned E-field','Single ended probe (P4)']
            legends = [['x','y','z'],'P4']
            ylabels = ['mV/m','V']
            
            if plot_options.style == '-':
                [axes[0].plot(Epoch_chunk[Mask_chunk[i] >= mask_limit],EF_in_xyz_calibr_chunk[i][Mask_chunk[i] >= mask_limit], lw = linewidth, label=legends[0][i]) for i in range(3)]
                axes[1].plot(Epoch_chunk[Mask_chunk[3] >= mask_limit],LP_diffs_chunk[3][Mask_chunk[3] >= mask_limit], lw = linewidth, label=legends[1])

            # elif plot_options.style == 'o':
            #     [axes[0].plot(Epoch_chunk[Mask_chunk[i] >= mask_limit],LP_diffs_chunk[i][Mask_chunk[i] >= mask_limit],'o', markersize = 1, linestyle = 'None',legends=labels[i]) for i in range(3)]
            #     axes[2].plot(Epoch_chunk[Mask_chunk[3] >= mask_limit],LP_diffs_chunk[3][Mask_chunk[3] >= mask_limit],'o', markersize = 1, linestyle = 'None',legends=labels[3])


            [axes[i].xaxis.set_major_formatter(FuncFormatter(dynamic_time_formatter)) for i in range(len(axes))]
            [axes[i].set_ylabel(ylabels[i]) for i in range(len(axes))]
            [axes[i].legend(loc='upper right') for i in range(len(axes))]
            [axes[i].set_title(titles[i]) for i in range(len(axes))]
            [axes[i].grid(True) for i in range(len(axes))]
            
            fig.canvas.manager.set_window_title("EF_xyz " + window_title)

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
    scatter_1 = px.scatter(x=Epoch, y=LP_diffs[0], labels={'y': 'P12'})
    scatter_2 = px.scatter(x=Epoch, y=LP_diffs[1], labels={'y': 'P23'})
    # scatter_3 = px.scatter(x=Epoch[lower_i:upper_i], y=LP_diffs[2][lower_i:upper_i], labels={'y': 'P34'})
    # scatter_4 = px.scatter(x=Epoch[lower_i:upper_i], y=LP_diffs[3][lower_i:upper_i], labels={'y': 'P4'})

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

