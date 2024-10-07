
# RPWI Electric field data
# By: Arvid Wennerström"


# ========================================================================
# ========================================================================
#                           IMPORTS
# ========================================================================
# ========================================================================
import numpy as np, scipy
from math import ceil

from rpwi_data import rpwi_data, coeffs_TM2voltage  
from plotting_support import datenum_to_tt2000, tt2000_to_readable, dynamic_time_formatter
# ========================================================================



# ========================================================================
# ========================================================================
#               INPUT OPTIONS AND SETTINGS FOR PLOTTING
# ========================================================================
# ========================================================================

#  Input overall directory path, where folders "datasets", "data_created" and "spice" exist
# ========================================================================
rootDir = "C:/Users/1/Onedrive - KTH/MEX/IRF" # My desktop
rootDir = "C:/Users/arvwe/Onedrive - KTH/MEX/IRF" # My laptop
rootDir = "C:/Users/arvidwen/Onedrive - KTH/MEX/IRF" # KTH computers



# Input date and time (optional) of the data to plot.
# ========================================================================#
date = '240820'

# Specify time period on format "HH:MM:SS.ns*"
# Precision is also optional, "08" and "08.30:20.500" are equally valid 
plot_time_period = ['21:07:53.629000','23:59']



# Select masking option 
# ========================================================================
# 1: TM saturation (either high or low).
# 2: Interference with RIME instrument.
# 3: Lunar wake data, will overwrite '1'.
# 4: High quality data, this does not contain any contaminations.
mask_limit = 4



# Select what data to plot and style
# ========================================================================
# NOTE: LP_data is now always converted to voltage from TM, so unit selection 
# does nothing. This might be desirable to change in the future.
plot_looks = {
    'Style':'lines',
    'Unit':'V'}

plots_2_show = {
    'Potential at probes': True,
    'Probe diffs': False,
    'E-field': False,
    'E-field_xyz': False,
    'DC_offset': True,
    'Sweep': False,
    'Plasmasphere idea': False}



# Select type of plots
# ========================================================================
plot_matplotlib = True 
plot_plotly = False



# Additional options, these will probably be removed in the future
# ========================================================================
DC_offset_option = 3


# Limit size of chunks to ensure better plot performance.
# I got 1024*512 (=524 288) from Ilona as a benchmark size, 
# but found that 5e6 works well for me
chunk_size_limit = 5e6  



# ========================================================================
# ========================================================================
#                           READ AND HANDLE DATA
# ========================================================================
# ========================================================================

#   Add SPICE-kernels, for converting time between tt2000 and human-readable
# ========================================================================
import spiceypy as spice    
ls_spice_kernels = [
    (rootDir + '/spice/JUICE/kernels/lsk/naif0012.tls'),
    (rootDir + '/spice/JUICE/kernels/sclk/juice_step_240828_v01.tsc')]
spice.furnsh(ls_spice_kernels)



#               Load E-field data
# ========================================================================
loaded_LP = np.load(rootDir + "/data_created/LP-SID1_" + date + ".npz")
Epoch = loaded_LP["Epoch"]
LP_diffs_TM = loaded_LP["LP_data"]
Mask = loaded_LP["Mask"]

date_long = '20' + date[0:2] + '-' + date[2:4] + '-' + date[4:6]



#               Load sweep data 
# NOTE: This should use regular expr.
# ========================================================================
if date == '240820':
    mat_data = scipy.io.loadmat(rootDir + '/sweep_data/lpstruc_20' + date + '_2_20241002.mat')
    t_sweep = mat_data['lpstruc']['t_sweep'][0][0]
    Usc = mat_data['lpstruc']['Usc'][0][0]
    Ni = abs(mat_data['lpstruc']['Ni'][0][0])

    Epoch_sweep = [datenum_to_tt2000(t[0])*1e9 for t in t_sweep]  


# Select only data of specific times, when a time interval has been provided
# ========================================================================
if plot_time_period != []:
    utc_time_begin = date_long + 'T' + plot_time_period[0]
    utc_time_end = date_long + 'T' + plot_time_period[1]

    # Convert UTC to Ephemeris Time (ET), ET to TT (Terrestrial Time)
    # and then multiply by 1e9 to get on same format as "Epoch"
    tt_time_begin = spice.unitim(spice.utc2et(utc_time_begin), "ET", "TT")*1e9
    tt_time_end = spice.unitim(spice.utc2et(utc_time_end), "ET", "TT")*1e9

    included_data = (Epoch > tt_time_begin) & (Epoch < tt_time_end)
    Epoch = Epoch[included_data]
    LP_diffs_TM = LP_diffs_TM[:,included_data]
    Mask = Mask[:,included_data]



# Conert raw data to voltages, from TM units
# ========================================================================
LP_diffs = coeffs_TM2voltage*LP_diffs_TM



# Get SC potential at each probe, by adding differentials.
# Use this to approximate DC offset coefficients and get
# calibrated potentials.
# ================================================
LP_potentials = np.zeros(LP_diffs.shape)
LP_potentials[3] = LP_diffs[3]
LP_potentials[2] = LP_potentials[3] + LP_diffs[2]
LP_potentials[1] = LP_potentials[2] + LP_diffs[1]
LP_potentials[0] = LP_potentials[1] + LP_diffs[0]



# Calculating linear DC offsets, using U4 = Un' = Un + dUn, for n = 1,2,3
# ================================================
# DC_offsets are np.arrays, with each value corresponding to a data value, to allow
# changes in offsets based on for example changes in plasma density. They are initialized below.
DC_offset1 = np.zeros(len(LP_potentials[3]))
DC_offset2 = np.zeros(len(LP_potentials[3]))
DC_offset3 = np.zeros(len(LP_potentials[3]))


# Very simple, just average over all of the data. Works poorly when data includes
# measurements in both high and low density plasma, such as plasmasphere on 20240820
if DC_offset_option == 1:
    DC_offset1 = DC_offset1 + np.mean(LP_potentials[3] - LP_potentials[0])
    DC_offset2 = DC_offset2 + np.mean(LP_potentials[3] - LP_potentials[1])
    DC_offset3 = DC_offset3 + np.mean(LP_potentials[3] - LP_potentials[2])


# A bit less simple: The offsets are averaged over sections of the data, 
# based on plot shape (voltage value corresponds to plasma density, 
# affecting the need for offset correction).
if DC_offset_option == 2:
    weight = 0.8
    # First chunk will be from start until first split, then  
    tt_splits = [
        spice.unitim(spice.utc2et(date_long), "ET", "TT")*1e9,
        spice.unitim(spice.utc2et(date_long + 'T20:38:10'), "ET", "TT")*1e9,
        spice.unitim(spice.utc2et(date_long + 'T20:13:10'), "ET", "TT")*1e9,
        spice.unitim(spice.utc2et(date_long + 'T21:30:30'), "ET", "TT")*1e9,
        spice.unitim(spice.utc2et(date_long + 'T22:33:40'), "ET", "TT")*1e9,
        spice.unitim(spice.utc2et(date_long + 'T23:03:40'), "ET", "TT")*1e9,
    ]
    for split_point in tt_splits:
        split = Epoch > split_point

        if split_point == spice.unitim(spice.utc2et(date_long + 'T21:30:30'), "ET", "TT")*1e9:        
            DC_offset1[split] = 0.1
            DC_offset2[split] = 0.1
            DC_offset3[split] = 0.05
        else:
            DC_offset1[split] = weight*np.mean(LP_potentials[3][split] - LP_potentials[0][split])
            DC_offset2[split] = weight*np.mean(LP_potentials[3][split] - LP_potentials[1][split])
            DC_offset3[split] = weight*np.mean(LP_potentials[3][split] - LP_potentials[2][split])



# Trying to relate DC offset to plasma density. This works pretty well for plasmasphere data, 
# but I don't know how well it relates to physics
if DC_offset_option == 3:
    f_Ni = scipy.interpolate.interp1d(Epoch_sweep,[i[0] for i in Ni], kind='linear', fill_value='extrapolate')
    N_interpolated = f_Ni(Epoch)

    N_max = round(max(Ni)[0])
    power = 0.2
    
    low_n_offset = [1, 0.8, 0.7]

    DC_offset1 = low_n_offset[0] + (-low_n_offset[0])/(N_max**power)*N_interpolated**power
    DC_offset2 = low_n_offset[1] + (-low_n_offset[1])/(N_max**power)*N_interpolated**power
    DC_offset3 = low_n_offset[2] + (-low_n_offset[2])/(N_max**power)*N_interpolated**power



# Apply DC offsets
# ================================================
LP_potentials_calibr = np.zeros(LP_potentials.shape)
LP_potentials_calibr[0] = LP_potentials[0] + DC_offset1
LP_potentials_calibr[1] = LP_potentials[1] + DC_offset2
LP_potentials_calibr[2] = LP_potentials[2] + DC_offset3
LP_potentials_calibr[3] = LP_potentials[3]



# Calculate new differentials using calibrated potentials
# ================================================
LP_diffs_calibr = np.zeros(LP_potentials_calibr.shape)
LP_diffs_calibr[3] = LP_potentials_calibr[3]
LP_diffs_calibr[2] = LP_potentials_calibr[2] - LP_potentials_calibr[3]
LP_diffs_calibr[1] = LP_potentials_calibr[1] - LP_potentials_calibr[2]
LP_diffs_calibr[0] = LP_potentials_calibr[0] - LP_potentials_calibr[1]



# Get E-field between probes
# ================================================
if plots_2_show['E-field']: 
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
if plots_2_show['E-field_xyz']:

    # M is transformation matrix in: U = M*E, where U is voltage and E is electric field
    M = [rpwi_data.LP12_distance,
         rpwi_data.LP23_distance,
         rpwi_data.LP34_distance]
    M_inv = np.linalg.inv(M)
    EF_in_xyz = 1e3*np.matmul(M_inv,LP_diffs[:3])
    EF_in_xyz_calibr = 1e3*np.matmul(M_inv,LP_diffs_calibr[:3])



# ================================================
# ================================================
#               Plot (using matplotlib)
# ================================================
# ================================================
if plot_matplotlib:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, MaxNLocator


    # Divide the data into chunks, by either number of data points in "limit" OR
    # if a time interval has been specified. Then plot a new figure for each chunk.
    # ================================================================
    if plot_time_period == []:
        number_of_data_chunks = ceil(len(Epoch)/chunk_size_limit)
    else:
        number_of_data_chunks = 1
    chunk_size = round(len(Epoch)/number_of_data_chunks)


    for chunk_i in range(number_of_data_chunks):
        # Decide and get data that belongs to this chunk, 
        # by creating a True/False masking array "chunk",
        # that tells which values to include.  
        # ================================================  
        lower_i = int(chunk_size*chunk_i)
        upper_i = int(chunk_size*(chunk_i+1)-1)
        chunk = np.zeros(Epoch.shape, dtype=bool)
        chunk[lower_i:upper_i] = True



        # General plot setup
        # ================================================
        window_title = date_long + ' Chunk ' + str(chunk_i+1) + ' (of ' + str(number_of_data_chunks) + ')'
        plot_style_options = [
            {
                'linewidth': '0.5',},
            
            {
                'marker': 'o',
                'markersize': '0.5',
                'linestyle': 'None'}]

        if plot_looks['Style'] == 'lines':
            plot_style = plot_style_options[0]
        elif plot_looks['Style'] == 'dots':
            plot_style = plot_style_options[1]
        xlabel = "Timestamp (UTC)"



        # Plot SC potential at each probe
        # ================================================
        if plots_2_show['Potential at probes']:
            legends = [
                ['P01 (P12+P23+P34+P04)', 'P02 (P23+P34+P04)', 'P03 (P34+P04)', 'P04'],
                ['U1 + dU1', 'U2 + dU2', 'U3 + dU3', 'U4']]  
            titles = [
                'Probe to SC potential',
                'Calibrated potentials']


            fig, axes = plt.subplots(2,1,sharex=True)
            [axes[0].plot(Epoch[chunk][Mask[i,chunk] >= mask_limit],LP_potentials[i,chunk][Mask[i,chunk] >= mask_limit], **plot_style, label=legends[0][i]) for i in range(4)]
            [axes[1].plot(Epoch[chunk][Mask[i,chunk] >= mask_limit],LP_potentials_calibr[i,chunk][Mask[i,chunk] >= mask_limit], **plot_style, label=legends[1][i]) for i in range(4)]
            

            # Setup for the subplots
            for i in range(len(axes)):
                axes[i].xaxis.set_major_formatter(FuncFormatter(dynamic_time_formatter))
                axes[i].set_xlabel(xlabel)
                axes[i].set_ylabel(plot_looks['Unit'])
                axes[i].set_title(titles[i])
                axes[i].legend(loc='upper right')
                axes[i].grid(True)  
            fig.canvas.manager.set_window_title("LP potentials " + window_title)


            # Rotate the x-axis labels for better readability
            plt.gcf().autofmt_xdate()
    


        # Plot probe differentials (either in voltages or TM)
        # ================================================
        if plots_2_show['Probe diffs']:
            legends = ['P12', 'P23', 'P34', 'P4']
            titles = [
                'Probe differentials',
                'Calibrated probe diff.',
                'Single ended (P4)']


            fig, axes = plt.subplots(3,1,sharex=True)
            [axes[0].plot(Epoch[chunk][Mask[i,chunk] >= mask_limit],LP_diffs[i,chunk][Mask[i,chunk] >= mask_limit], **plot_style, label=legends[i]) for i in range(3)]
            [axes[1].plot(Epoch[chunk][Mask[i,chunk] >= mask_limit],LP_diffs_calibr[i,chunk][Mask[i,chunk] >= mask_limit], **plot_style, label=legends[i]) for i in range(3)]          
            axes[2].plot(Epoch[chunk][Mask[3,chunk] >= mask_limit],LP_diffs_calibr[3,chunk][Mask[3,chunk] >= mask_limit], **plot_style, label=legends[3])


            # Setup for the subplots
            for i in range(len(axes)):
                axes[i].xaxis.set_major_formatter(FuncFormatter(dynamic_time_formatter))
                axes[i].set_xlabel(xlabel)
                axes[i].set_ylabel(plot_looks['Unit'])
                axes[i].set_title(titles[i])
                axes[i].legend(loc='upper right')
                axes[i].grid(True)  
            fig.canvas.manager.set_window_title('LP differentials ' + window_title)

            # Rotate the x-axis labels for better readability
            plt.gcf().autofmt_xdate()



        # Plot electric field strength in probe differentials
        # ================================================
        if plots_2_show['E-field']:
            fig, axes = plt.subplots(3,1, sharex=True)
            titles = [
                'Probe differential E-field',
                'Calibrated diff. E-field',
                'Single ended probe (P4)']
            legends = ['P12', 'P23', 'P34','P4']
            ylabels = ['mV/m','mV/m','V']


            [axes[0].plot(Epoch[chunk][Mask[i,chunk] >= mask_limit],EF_in_LPs[i,chunk][Mask[i,chunk] >= mask_limit], **plot_style, label=legends[i]) for i in range(3)]
            [axes[1].plot(Epoch[chunk][Mask[i,chunk] >= mask_limit],EF_in_LPs_calibr[i,chunk][Mask[i,chunk] >= mask_limit], **plot_style, label=legends[i]) for i in range(3)]
            axes[2].plot(Epoch[chunk][Mask[3,chunk] >= mask_limit],LP_diffs[3,chunk][Mask[3,chunk] >= mask_limit], **plot_style, label=legends[3])

            for i in range(len(axes)):
                axes[i].xaxis.set_major_formatter(FuncFormatter(dynamic_time_formatter))
                axes[i].set_ylabel(ylabels[i])
                axes[i].legend(loc='upper right')
                axes[i].set_title(titles[i])
                axes[i].grid(True)

            fig.canvas.manager.set_window_title("EF " + window_title)

            # Rotate the x-axis labels for better readability
            plt.gcf().autofmt_xdate()



        # Plot electric field strenth in x,y,z-axes
        # ================================================
        if plots_2_show['E-field_xyz']:
            fig, axes = plt.subplots(2,1, sharex=True)
            titles = ['Calibrated axis aligned E-field','Single ended probe (P4)']
            legends = [['x','y','z'],'P4']
            ylabels = ['mV/m','V']
            

            [axes[0].plot(Epoch[chunk][Mask[i,chunk] >= mask_limit],EF_in_xyz_calibr[i,chunk][Mask[i,chunk] >= mask_limit], **plot_style, label=legends[0][i]) for i in range(3)]
            axes[1].plot(Epoch[chunk][Mask[3,chunk] >= mask_limit],LP_diffs[3,chunk][Mask[3,chunk] >= mask_limit], **plot_style, label=legends[1])

            for i in range(len(axes)):
                axes[i].xaxis.set_major_formatter(FuncFormatter(dynamic_time_formatter))
                axes[i].set_ylabel(ylabels[i])
                axes[i].legend(loc='upper right')
                axes[i].set_title(titles[i])
                axes[i].grid(True)
            
            fig.canvas.manager.set_window_title("EF_xyz " + window_title)

            # Rotate the x-axis labels for better readability
            plt.gcf().autofmt_xdate()
            


        # DC offset used in calibration
        # ========================================================================
        if plots_2_show['DC_offset']:
            plt.figure()
            plt.plot(Epoch,DC_offset1, label='DC_offset1')
            plt.plot(Epoch,DC_offset2, label='DC_offset2')
            plt.plot(Epoch,DC_offset3, label='DC_offset3')
            plt.legend(loc='upper right')
            plt.grid(True)
            plt.gca().xaxis.set_major_formatter(FuncFormatter(dynamic_time_formatter))
            plt.gcf().autofmt_xdate()




        # Plot sweep SC potential
        # ========================================================================
        if plots_2_show['Sweep']:    

            # Interpolation of low-resolution data to high-resolution time points
            f = scipy.interpolate.interp1d(Epoch_sweep,[u[0] for u in Usc], kind='linear', fill_value='extrapolate')
            Usc_interpolated = f(Epoch)
    

            fig, axes = plt.subplots(2,1, sharex=True)
            titles = [
                'E-field and sweep SC potentials',
                'Difference']

            axes[0].plot(Epoch,LP_potentials[3], **plot_style, label='P04')
            axes[0].plot(Epoch,-Usc_interpolated, **plot_style, label='-Usc')
            axes[1].plot(Epoch,LP_potentials[3]+Usc_interpolated, **plot_style, label='P04 + Usc')

            for i in range(len(axes)):
                axes[i].xaxis.set_major_locator(MaxNLocator(nbins=10))
                axes[i].xaxis.set_major_formatter(FuncFormatter(dynamic_time_formatter))
                axes[i].set_ylabel('V')
                axes[i].legend(loc='upper right')
                axes[i].set_title(titles[i])
                axes[i].grid(True)
            fig.canvas.manager.set_window_title("Sweep SC potential " + window_title)

            # Rotate the x-axis labels for better readability
            plt.gcf().autofmt_xdate()



        # ========================================================================
        if plots_2_show['Plasmasphere idea']:
            U1 = LP_potentials[0,chunk][Mask[0,chunk] >= mask_limit]
            U2 = LP_potentials[1,chunk][Mask[1,chunk] >= mask_limit]
            U3 = LP_potentials[2,chunk][Mask[2,chunk] >= mask_limit]
            U4 = LP_potentials[3,chunk][Mask[3,chunk] >= mask_limit]
            
            fig, axes = plt.subplots(2,1,sharex=True)
            axes[0].plot(U1,U4, **plot_style_options[1], label='U1')
            axes[0].plot(U2,U4, **plot_style_options[1], label='U2')
            axes[0].plot(U3,U4, **plot_style_options[1], label='U3')
            axes[1].plot(U4-U1,U4, **plot_style_options[1], label='U4-U1')
            axes[1].plot(U4-U2,U4, **plot_style_options[1], label='U4-U2')
            axes[1].plot(U4-U3,U4, **plot_style_options[1], label='U4-U3')

            for i in range(2):
                axes[i].legend(loc='lower right')
                axes[i].set_ylabel('U4')   
                axes[i].grid(True)
                axes[i].set_title('Alternative ' + str(i))
            fig.canvas.manager.set_window_title("Plasmasphere idea " + window_title)


        # This chunk's plots
        # ================================================
        plt.show()





#                 Plot (using plotly)
# ================================================
if plot_plotly:
    pass
   
    # import plotly.graph_objects as go
    # import plotly.express as px
    # from plotly.subplots import make_subplots

    # fig = make_subplots(rows=1, cols=1, shared_xaxes=True)

    # Create plotly express scatter plots
    # scatter_1 = px.scatter(x=Epoch, y=LP_diffs[0], labels={'y': 'P12'})
    # scatter_2 = px.scatter(x=Epoch, y=LP_diffs[1], labels={'y': 'P23'})
    # scatter_3 = px.scatter(x=Epoch[lower_i:upper_i], y=LP_diffs[2][lower_i:upper_i], labels={'y': 'P34'})
    # scatter_4 = px.scatter(x=Epoch[lower_i:upper_i], y=LP_diffs[3][lower_i:upper_i], labels={'y': 'P4'})

    # Add scatter plots to subplots
    # fig.add_trace(scatter_1.data[0], row=1, col=1)
    # fig.add_trace(scatter_2.data[0], row=1, col=1)
    # fig.add_trace(scatter_3.data[0], row=1, col=1)
    # fig.add_trace(scatter_4.data[0], row=2, col=1)

    # Update layout for the figure
    # fig.update_layout(title_text="LP Data with Shared X-Axis")

    # Show grid for both subplots
    # fig.update_xaxes(showgrid=True)
    # fig.update_yaxes(showgrid=True)

    # Add x-axis and y-axis labels
    # fig.update_xaxes(title_text="Time", row=2, col=1)
    # fig.update_yaxes(title_text=ylabel[0], row=1, col=1)
    # fig.update_yaxes(title_text=ylabel[1], row=2, col=1)

    # Show the plot
    # fig.show()



    # # Melt the DataFrame to a long format
    # df_melted = df.melt(id_vars='x', value_vars=['y1', 'y2', 'y3', 'y4'], var_name='Series', value_name='y')


    # # Create scatter plot
    # fig = px.scatter(df_melted, x='x', y ='y', color='Series')

    # # Show the plot
    # fig.show()



# # NOTE: This has nothing to do with plotting, just a request by Jan-Erik
# # ================================================
# if date == '240820':
#     import scipy.io
#     Epoch_test = Epoch[chunk][Mask[3,chunk] >= mask_limit]
#     P4_test = LP_potentials[3,chunk][Mask[3,chunk] >= mask_limit]
#     scipy.io.savemat((rootDir + '/P04_20240820T200000.mat'),{'Epoch': Epoch_test,'P04': P4_test})

#     fig, ax = plt.subplots(1,1)
#     ax.plot(Epoch_test,P4_test)
#     ax.xaxis.set_major_formatter(FuncFormatter(dynamic_time_formatter))
#     plt.show()