
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
from plotting_support import Data_struct,datenum_to_tt2000, dynamic_time_formatter
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
date = '240823'
date_long = '20' + date[0:2] + '-' + date[2:4] + '-' + date[4:6]


# Specify time period on format "HH:MM:SS.ns*"
# Precision is also optional, "08" and "08.30:20.500" are equally valid 
plot_time_period = []

# plot_time_period = ['21:07:53.629000','23:59'] # This is the period for which we have sweep data on 240820
plot_time_period = ['04','05:30']


# Select masking option 
# ========================================================================
# 1: TM saturation (either high or low).
# 2: Interference with RIME instrument.
# 3: Lunar wake data, will overwrite '1'.
# 4: High quality data, this does not contain any contaminations.
mask_limit = 4



# Additional options, these will probably be removed in the future
# ========================================================================
DC_offset_option = 1


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



# Load E-field data and select only data of specified time, if any is provided
# ========================================================================
loaded_LP = np.load(rootDir + "/data_created/LP-SID1_" + date + ".npz")
Epoch = loaded_LP["Epoch"]
Mask = loaded_LP["Mask"]
LP_data = loaded_LP["LP_data"]


# If a time period has been specified, only extract the data of said period. 
if plot_time_period != []:
    utc_time_begin = date_long + 'T' + plot_time_period[0]
    utc_time_end = date_long + 'T' + plot_time_period[1]

    # Convert UTC to Ephemeris Time (ET), ET to TT (Terrestrial Time)
    # and then multiply by 1e9 to get on same format as "Epoch"
    tt_time_begin = spice.unitim(spice.utc2et(utc_time_begin), "ET", "TT")*1e9
    tt_time_end = spice.unitim(spice.utc2et(utc_time_end), "ET", "TT")*1e9

    included_data = (Epoch > tt_time_begin) & (Epoch < tt_time_end)
    Epoch = Epoch[included_data]
    LP_data = LP_data[:,included_data]
    Mask = Mask[:,included_data]


# All data is saved in these "Data_struct"-objects, found in "plotting_support.py".
# Among the actual data, they contain epoch time, masking and units etc. complementing the data. 
LP_diffs_TM = Data_struct(LP_data[0:3],Epoch,Mask,'TM',['P12', 'P23', 'P34'],'Probe differentials')
Single_ended_TM = Data_struct([LP_data[3]],Epoch,Mask,'TM',['P4'],'Single ended (P04) potential')


#               Load sweep data 
# NOTE: This should use regular expr.
# ========================================================================
if date == '240820':
    mat_data = scipy.io.loadmat(rootDir + '/sweep_data/lpstruc_20' + date + '_2_20241002.mat')
    t_sweep = mat_data['lpstruc']['t_sweep'][0][0]
    Usc = mat_data['lpstruc']['Usc'][0][0]
    Ni = abs(mat_data['lpstruc']['Ni'][0][0])

    Epoch_sweep = [datenum_to_tt2000(t[0])*1e9 for t in t_sweep]  


# Convert raw data to voltages, from TM units
# ========================================================================
LP_diffs = Data_struct(coeffs_TM2voltage[0:3]*LP_diffs_TM.data,Epoch,Mask,'V',['P12', 'P23', 'P34'],'Probe differentials')
Single_ended = Data_struct(Single_ended_TM.data*coeffs_TM2voltage[3],Epoch,Mask,'V',['P4'],'Single ended (P4) potential')


# Get SC potential at each probe, by adding differentials.
# Use this to approximate DC offset coefficients and get
# calibrated potentials.
# ================================================
LP_potentials = Data_struct(np.zeros((4,len(Epoch))),Epoch,Mask,'V',['P01 (P12+P23+P34+P04)', 'P02 (P23+P34+P04)', 'P03 (P34+P04)', 'P04'],'Probe potentials')
LP_potentials.data[3] = Single_ended.data
LP_potentials.data[2] = LP_potentials.data[3] + LP_diffs.data[2]
LP_potentials.data[1] = LP_potentials.data[2] + LP_diffs.data[1]
LP_potentials.data[0] = LP_potentials.data[1] + LP_diffs.data[0]



# Calculating linear DC offsets, using U4 = Un' = Un + dUn, for n = 1,2,3
# ================================================
# DC_offsets are np.arrays, with each value corresponding to a data value, to allow
# changes in offsets based on for example changes in plasma density. They are initialized below.
DC_offset = Data_struct(np.zeros((3,len(Epoch))),Epoch,Mask,'V',['Offset 1','Offset 2','Offset 3'],'DC offsets')

DC_offset1 = np.zeros(len(LP_potentials.data[3]))
DC_offset2 = np.zeros(len(LP_potentials.data[3]))
DC_offset3 = np.zeros(len(LP_potentials.data[3]))


# Very simple, just average over all of the data. Works poorly when data includes
# measurements in both high and low density plasma, such as plasmasphere on 20240820
if DC_offset_option == 1:
    DC_offset.data[0] = DC_offset.data[0] + np.mean(LP_potentials.data[3] - LP_potentials.data[0])
    DC_offset.data[1] = DC_offset.data[1] + np.mean(LP_potentials.data[3] - LP_potentials.data[1])
    DC_offset.data[2] = DC_offset.data[2] + np.mean(LP_potentials.data[3] - LP_potentials.data[2])

    # DC_offset1 = DC_offset1 + np.mean(LP_potentials.data[3] - LP_potentials.data[0])
    # DC_offset2 = DC_offset2 + np.mean(LP_potentials.data[3] - LP_potentials.data[1])
    # DC_offset3 = DC_offset3 + np.mean(LP_potentials.data[3] - LP_potentials.data[2])


# # A bit less simple: The offsets are averaged over sections of the data, 
# # based on plot shape (voltage value corresponds to plasma density, 
# # affecting the need for offset correction).
# if DC_offset_option == 2:
#     weight = 0.8
#     # First chunk will be from start until first split, then  
#     tt_splits = [
#         spice.unitim(spice.utc2et(date_long), "ET", "TT")*1e9,
#         spice.unitim(spice.utc2et(date_long + 'T20:38:10'), "ET", "TT")*1e9,
#         spice.unitim(spice.utc2et(date_long + 'T20:13:10'), "ET", "TT")*1e9,
#         spice.unitim(spice.utc2et(date_long + 'T21:30:30'), "ET", "TT")*1e9,
#         spice.unitim(spice.utc2et(date_long + 'T22:33:40'), "ET", "TT")*1e9,
#         spice.unitim(spice.utc2et(date_long + 'T23:03:40'), "ET", "TT")*1e9,
#     ]
#     for split_point in tt_splits:
#         split = Epoch > split_point

#         if split_point == spice.unitim(spice.utc2et(date_long + 'T21:30:30'), "ET", "TT")*1e9:        
#             DC_offset1[split] = 0.1
#             DC_offset2[split] = 0.1
#             DC_offset3[split] = 0.05
#         else:
#             DC_offset1[split] = weight*np.mean(LP_potentials.data[3][split] - LP_potentials.data[0][split])
#             DC_offset2[split] = weight*np.mean(LP_potentials.data[3][split] - LP_potentials.data[1][split])
#             DC_offset3[split] = weight*np.mean(LP_potentials.data[3][split] - LP_potentials.data[2][split])


# # Trying to relate DC offset to plasma density. This works pretty well for plasmasphere data, 
# # but I don't know how well it relates to physics
# if DC_offset_option == 3:
#     f_Ni = scipy.interpolate.interp1d(Epoch_sweep,[i[0] for i in Ni], kind='linear', fill_value='extrapolate')
#     N_interpolated = f_Ni(Epoch)

#     N_max = max(Ni)[0]
#     power = 1
    
#     low_n_offset = [1, 0.8, 0.7]

#     DC_offset1 = low_n_offset[0] + (-low_n_offset[0])/(N_max**power)*N_interpolated**power
#     DC_offset2 = low_n_offset[1] + (-low_n_offset[1])/(N_max**power)*N_interpolated**power
#     DC_offset3 = low_n_offset[2] + (-low_n_offset[2])/(N_max**power)*N_interpolated**power


# # Trying to relate DC offset to plasma density. This works pretty well for plasmasphere data, 
# # but I don't know how well it relates to physics
# if DC_offset_option == 4:
#     low_n_offset = [1, 0.8, 0.7]
#     U4_max = max(LP_potentials[3])

#     DC_offset1 = low_n_offset[0] + (-low_n_offset[0])/(U4_max)*LP_potentials[3]
#     DC_offset2 = low_n_offset[1] + (-low_n_offset[1])/(U4_max)*LP_potentials[3]
#     DC_offset3 = low_n_offset[2] + (-low_n_offset[2])/(U4_max)*LP_potentials[3]



# Apply DC offsets
# ================================================
LP_potentials_offset = Data_struct(np.zeros((4,len(Epoch))),Epoch,Mask,'V',['U1 + dU1', 'U2 + dU2', 'U3 + dU3', 'U4'],'Offset probe potentials')
LP_potentials_offset.data[3] = LP_potentials.data[3]
LP_potentials_offset.data[2] = LP_potentials.data[2] + DC_offset.data[2]
LP_potentials_offset.data[1] = LP_potentials.data[1] + DC_offset.data[1]
LP_potentials_offset.data[0] = LP_potentials.data[0] + DC_offset.data[0]



# Calculate new differentials using calibrated potentials
# ================================================
LP_diffs_offset = Data_struct(np.zeros((3,len(Epoch))),Epoch,Mask,'V',['P12','P23','P34'],'Offset probe differnetials')
LP_diffs_offset.data[2] = LP_potentials_offset.data[2] - Single_ended.data
LP_diffs_offset.data[1] = LP_potentials_offset.data[1] - LP_potentials_offset.data[2]
LP_diffs_offset.data[0] = LP_potentials_offset.data[0] - LP_potentials_offset.data[1]



# Get E-field between probes
# ================================================ 
# 1e3's are used to give EF in mV/m, rather than V/m
# NOTE: This is for mode 0
EF_in_LPs = Data_struct(np.zeros((3,len(Epoch))),Epoch,Mask,'mV/m',['P12', 'P23', 'P34'],'Electric field')
EF_in_LPs.data[0] = 1e3*LP_diffs.data[0]/np.linalg.norm(rpwi_data.LP12_distance)
EF_in_LPs.data[1] = 1e3*LP_diffs.data[1]/np.linalg.norm(rpwi_data.LP23_distance)
EF_in_LPs.data[2] = 1e3*LP_diffs.data[2]/np.linalg.norm(rpwi_data.LP34_distance)

EF_in_LPs_offset = Data_struct(np.zeros((3,len(Epoch))),Epoch,Mask,'mV/m',['P12', 'P23', 'P34'],'Electric field offset')
EF_in_LPs_offset.data[0] = 1e3*LP_diffs_offset.data[0]/np.linalg.norm(rpwi_data.LP12_distance)
EF_in_LPs_offset.data[1] = 1e3*LP_diffs_offset.data[1]/np.linalg.norm(rpwi_data.LP23_distance)
EF_in_LPs_offset.data[2] = 1e3*LP_diffs_offset.data[2]/np.linalg.norm(rpwi_data.LP34_distance)



# Calculate E-field in SC x-,y- and z-axes
# ================================================
# M is transformation matrix in: U = M*E, where U is voltage and E is electric field
M = [rpwi_data.LP12_distance,
    rpwi_data.LP23_distance,
    rpwi_data.LP34_distance]
M_inv = np.linalg.inv(M)
EF_in_xyz = Data_struct(1e3*np.matmul(M_inv,LP_diffs.data[:3]),Epoch,Mask,'mV/m',['x','y','z'],'EF in xyz')
EF_in_xyz_offset = Data_struct(1e3*np.matmul(M_inv,LP_diffs_offset.data[:3]),Epoch,Mask,'mV/m',['x','y','z'],'EF in xyz')



# ================================================
# ================================================
#               Plot (using matplotlib)
# ================================================
# ================================================
if True:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, MaxNLocator


    # Divide the data into chunks, by either number of data points in "chunk_size_limit" OR
    # if a time interval has been specified. Then plot a new figure for each chunk.
    # ================================================================
    if plot_time_period == []:
        number_of_data_chunks = ceil(len(Epoch)/chunk_size_limit)
    else:
        number_of_data_chunks = 1
 

    for chunk_i in range(number_of_data_chunks):
        # Decide and get data that belongs to this chunk, 
        # by creating a True/False masking array "chunk",
        # that tells which values to include. 
        # ================================================ 
        chunk = np.array_split(np.zeros(Epoch.shape, dtype=bool),number_of_data_chunks)
        chunk[chunk_i][:] = True
        chunk = np.concatenate(chunk)


        # General plot setup
        # ================================================
        if number_of_data_chunks == 1:
            window_title = date_long
        else:
            window_title = date_long + ' Chunk ' + str(chunk_i+1) + ' (of ' + str(number_of_data_chunks) + ')'
        
        plot_styles = {
            'lines':{
                'linewidth': '0.5',},
            'dots':{
                'marker': 'o',
                'markersize': '0.5',
                'linestyle': 'None'}}


        # # Plot SC potential at each probe
        # # ================================================
        # fig, axes = plt.subplots(2,1,sharex=True)
        # LP_potentials.plot(axes[0],chunk,plot_styles['lines'],4)
        # LP_potentials_offset.plot(axes[1],chunk,plot_styles['lines'],4)
        # fig.canvas.manager.set_window_title("LP potentials " + window_title)



        # # Plot probe differentials (either in voltages or TM)
        # # ================================================
        # fig, axes = plt.subplots(3,1,sharex=True)
        # LP_diffs.plot(axes[0],chunk,plot_styles['lines'],4)
        # LP_diffs_offset.plot(axes[1],chunk,plot_styles['lines'],4)
        # Single_ended.plot(axes[2],chunk,plot_styles['lines'],4)
        # fig.canvas.manager.set_window_title("LP differentials " + window_title)



        # # Plot electric field strength in probe differentials
        # # ================================================
        # fig, axes = plt.subplots(3,1,sharex=True)
        # EF_in_LPs.plot(axes[0],chunk,plot_styles['lines'],4)
        # EF_in_LPs_offset.plot(axes[1],chunk,plot_styles['lines'],4)
        # Single_ended.plot(axes[2],chunk,plot_styles['lines'],4)
        # fig.canvas.manager.set_window_title("EFs " + window_title)



        # # Plot electric field strenth in x,y,z-axes
        # # ================================================
        # fig, axes = plt.subplots(3,1,sharex=True)
        # EF_in_xyz.plot(axes[0],chunk,plot_styles['lines'],4)
        # EF_in_xyz_offset.plot(axes[1],chunk,plot_styles['lines'],4)
        # Single_ended.plot(axes[2],chunk,plot_styles['lines'],4)
        # fig.canvas.manager.set_window_title("EFs in XYZ " + window_title)



        # DC offset used in calibration
        # ========================================================================
        fig = plt.figure()
        DC_offset.plot(fig.gca(),chunk,plot_styles['lines'])
        


        if False:
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



            # Jan-Erik's plasmasphere idea
            # ========================================================================
            if plots_2_show['Plasmasphere idea'] and date == '240820':
                U1 = LP_potentials[0,chunk][Mask[0,chunk] >= mask_limit]
                U2 = LP_potentials[1,chunk][Mask[1,chunk] >= mask_limit]
                U3 = LP_potentials[2,chunk][Mask[2,chunk] >= mask_limit]
                U4 = LP_potentials[3,chunk][Mask[3,chunk] >= mask_limit]
                
                fig, axes = plt.subplots(2,1,sharex=True)
                axes[0].plot(U1,U4, **plot_styles['dots'], label='U1')
                axes[0].plot(U2,U4, **plot_styles['dots'], label='U2')
                axes[0].plot(U3,U4, **plot_styles['dots'], label='U3')
                axes[1].plot(U4-U1,U4, **plot_styles['dots'], label='U4-U1')
                axes[1].plot(U4-U2,U4, **plot_styles['dots'], label='U4-U2')
                axes[1].plot(U4-U3,U4, **plot_styles['dots'], label='U4-U3')

                for i in range(2):
                    axes[i].legend(loc='lower right')
                    axes[i].set_ylabel('U4')   
                    axes[i].grid(True)
                    axes[i].set_title('Alternative ' + str(i))
                fig.canvas.manager.set_window_title("Plasmasphere idea " + window_title)



            # Plasma density in plasmasphere 
            # ========================================================================
            if plots_2_show['Plasma density'] and date == '240820':
                plt.figure()
                plt.plot(Epoch_sweep,Ni,**plot_styles['lines'])
                plt.yscale('log')
                plt.title('Plasma density')
                plt.ylabel('Number density [cm-3]')
                plt.gca().xaxis.set_major_formatter(FuncFormatter(dynamic_time_formatter))
                plt.grid(True)
                plt.gcf().autofmt_xdate()

                print(max(Ni),min(Ni))

        # This chunk's plots
        # ================================================
        plt.show()

