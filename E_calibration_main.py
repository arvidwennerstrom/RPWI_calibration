"""
RPWI Electric field data
By: Arvid Wennerström


                    DOCUMENTATION
========================================================================
========================================================================

Most of the data is stored in "Struct"-objects (found in "plotting_support.py"). 
These contain: 
    - data, an numpy.array of optional size
    - time, an array of time stamps in tt2000, which maps to the data points
    - mask, equal shape as "data", giving each data point a mapping value
    - unit of the data
    and more...

The default naming convention is: 
    "*Data type*_*Source*_*Coordinates*

    
Data types include [default unit]:
    B       -   Magnetic field [nT]
    E       -   Electric field [mV/m]
    v       -   Juice's velocity [km/s]
    pos     -   Juice's position [R_E]

    
Sources include:
    LP      -   Measured by Juice LPs
    IGRF    -   IGRF magnetic field model
    JMAG    - 

    
Coordinates include:
    LPs     -   Along LP differentials axes.
    SC      -   Juice's x-, y- and z-axes.
    GSE     -   Geocentric Solar Ecliptic.
    GSM     -   Geocentric Solar Magnetospheric.


For example: 
    "E_LP_SC" contains E-field data, as observed by Juice, in Juice's x,y,z.

Exceptions to this are often self explanatory, such as: 
    LP_potentials - Potential of the 4 Langmuir probes

"""



"""     IMPORTS
========================================================================
========================================================================
"""
import datetime, numpy as np, scipy, matplotlib.pyplot as plt
from E_calibr_support import Struct, datenum_to_tt2000, rpwi_data, coeffs_TM2voltage  
from IGRF import IGRF_B_field
from fit_to_E_obs import fit_to_E_obs



"""     INPUT OPTIONS AND SETTINGS FOR PLOTTING
========================================================================           
========================================================================
"""
# Input overall directory path, where folders "datasets", 
# "data_created" and "spice" exist
# ========================================================================
rootDir = "C:/Users/1/Onedrive - KTH/MEX/IRF" # My desktop
# rootDir = "C:/Users/arvwe/Onedrive - KTH/MEX/IRF" # My laptop
rootDir = "C:/Users/arvidwen/Onedrive - KTH/MEX/IRF" # KTH computers


# Input date and time (optional) of the data to plot.
# ========================================================================
# Specify on format: "YYYY-MM-DDTHH:MM:SS.ns*". 
# Date is mandatory, adding T is optional: "2024-08-20" and 2024-08-20T10" are both valid. 
# Time precision is optional: "*T08" and "*T08.30:20.500" are both valid. 
time_period = ['2024-08-23T04', '2024-08-23T05:30']

start_date = datetime.datetime.strptime(time_period[0].split('T')[0], "%Y-%m-%d")
end_date = datetime.datetime.strptime(time_period[1].split('T')[0], "%Y-%m-%d")
date_list = [(start_date + datetime.timedelta(days=x)).strftime("%Y%m%d") for x in range((end_date - start_date).days + 1)]


plot_styles = {
    'lines':{
        'linewidth': '1',},
    'lines_thin':{
        'linewidth': '0.5',},
    'dots':{
        'marker': 'o',
        'markersize': '0.5',
        'linestyle': 'None'}}




"""     LOAD LP DATA
========================================================================                        
========================================================================
The data is read from .npz files created by "data_creation.py". These contain an
Epoch array of time stamps and the LP measurements values for these, in TM units. 
The LP data has been handled to remove noise etc., which is handled by the 
"Mask"-array:
    1: TM saturation (either high or low).
    2: Interference with RIME instrument.
    3: Lunar wake data, will overwrite '1'.
    4: High quality data, does not contain any contaminations.

"""

#   Add SPICE-kernels, for converting time between Epoch_LEGA and human-readable
# ========================================================================
import spiceypy as spice    
ls_spice_kernels = [
    (rootDir + '/spice/JUICE/kernels/lsk/naif0012.tls'),
    (rootDir + '/spice/JUICE/kernels/sclk/juice_step_240828_v01.tsc')]
spice.furnsh(ls_spice_kernels)



# Load E-field data from all days included in the specified time period
# ========================================================================
i = 1
for date in date_list:  
    # Ensures loading of .npz files only for days where data exists
    try:
        loaded_LP = np.load(rootDir + "/data_created/LP-SID1_" + date + ".npz")
        
        # Set up data structures
        if i == 1:
            Epoch = loaded_LP["Epoch"]
            LP_data = loaded_LP["LP_data"]
            Mask = loaded_LP["Mask"]

        # Add to data
        else:
            Epoch = np.concatenate((Epoch,loaded_LP["Epoch"]),axis = 0)
            LP_data = np.concatenate((LP_data,loaded_LP["LP_data"]),axis = 1)
            Mask = np.concatenate((Mask,loaded_LP["Mask"]),axis = 1)

        i += 1
    except:
        pass
    


# Truncate data to only include specified time period
# ========================================================================
# Convert UTC to Ephemeris Time (ET), ET to TT (Terrestrial Time)
# and then multiply by 1e9 to get on same format as "Epoch" (tt2000)
included_data_LP = (Epoch > spice.unitim(spice.utc2et(time_period[0]), "ET", "TT")*1e9) & (Epoch < spice.unitim(spice.utc2et(time_period[1]), "ET", "TT")*1e9)
Epoch = Epoch[included_data_LP]
LP_data = LP_data[:,included_data_LP]
Mask = Mask[:,included_data_LP]



# Save TM data 
# Convert raw data to voltages, from TM units
# ========================================================================
LP_diffs_TM = Struct(LP_data,Epoch,Mask,'TM',['P12', 'P23', 'P34','P4 (Single ended)'],'Probe differentials')
LP_diffs = Struct(coeffs_TM2voltage*LP_diffs_TM.data,Epoch,Mask,'V',['P12', 'P23', 'P34','P04 (Single ended)'],'Probe differentials SID1')
LP_diffs.data[3] = -LP_diffs[3]


# Get SC potential at each probe, by adding differentials.
# ================================================
LP_potentials = Struct(np.zeros((4,len(Epoch))),Epoch,Mask,'V',['P01 (P12+P23+P34+P04)', 'P02 (P23+P34+P04)', 'P03 (P34+P04)', 'P04'],'Probe potentials')
LP_potentials.data[3] = LP_diffs[3]
LP_potentials.data[2] = LP_potentials.data[3] + LP_diffs.data[2]
LP_potentials.data[1] = LP_potentials.data[2] + LP_diffs.data[1]
LP_potentials.data[0] = LP_potentials.data[1] + LP_diffs.data[0]



"""     MEASURED ELECTRIC FIELD
======================================================================== 
========================================================================
"""

# Get E-field between probes
# ================================================ 
# 1e3's are used to give EF in mV/m, rather than V/m
E_LP_LPs = Struct(np.zeros((3,len(Epoch))),Epoch,Mask,'mV/m',['P12', 'P23', 'P34'],'E-field in LP differentials')
for i in range(3):
    E_LP_LPs.data[i] = 1e3*LP_diffs[i]/np.linalg.norm(rpwi_data.M_E2U[i]) 
    

# Calculate E-field in SC x-,y- and z-axes
# ================================================
E_LP_SC = Struct(1e3*np.matmul(rpwi_data.M_U2E,LP_diffs.data[:3]),Epoch,Mask,'mV/m',['Ex','Ey','Ez'],"E-field from LP's in SC coords")




"""     Load B- and v-data
======================================================================== 
========================================================================
"""
# Load OMNI2 data (V_sw and B-field)
# ========================================================================
from omni import load_omni_data
B_OMNI_GSE, SW_OMNI =  load_omni_data(rootDir,Epoch)



# Load J-MAG data (B-field)
# ========================================================================
from jmag import load_jmag_data
B_JMAG_SC = load_jmag_data(rootDir,Epoch) 



# Load IGRF data (B-field, E-field etc.) 
# ========================================================================
B_IGRF_GSM, B_IGRF_SC, E_IGRF_GSM, E_IGRF_SC, juice_V_SC = IGRF_B_field(rootDir,Epoch)



# Calculate -v x B 
# ========================================================================

# Only take x,y,z-components (not magnitude)
B = B_JMAG_SC.data[0:3,:]
v = np.zeros_like(B); v[0] = SW_OMNI.data[0]


# Get E-field (in mV/m)
E_JMAG_SC = np.zeros_like(B_JMAG_SC.data)
E_JMAG_SC[0:3,:] = np.cross(-v.T,B.T).T*1e-3
E_JMAG_SC[3,:] = np.linalg.norm(E_JMAG_SC, axis=0)
E_JMAG_SC = Struct(E_JMAG_SC,Epoch,None,'mV/M',['Ex','Ey','Ez','|E|'],'E-field from J-MAG in SC coords')




"""     Calculate E-field calibration 
========================================================================
========================================================================
"""
cross_coefficients_E = False


# NOTE: Can improve this, for example by looking at time period
if False:
    # For plasmasphere period, where B comes from IGRF and v from juice velocity w.r.t the Earth.
    # ================================================
    E_cal_SC, cal_coeff_a, cal_coeff_b, cal_coeff_c, coeff_Epoc = fit_to_E_obs(Epoch,E_LP_SC.data,E_IGRF_SC.data,LP_potentials[3],cross_coefficients_E)
    E_cal_SC = Struct(E_cal_SC,Epoch,None,'mV/m',['x - cal.','y - cal.','z - cal.'],'Calibrated E-field, from IGRF')


if True:
    # For solar wind period, where B comes from J-MAG and v = v_sw comes from ACE (OMNI2). 
    # ================================================
    E_cal_SC, cal_coeff_a, cal_coeff_b, cal_coeff_c, coeff_Epoc = fit_to_E_obs(Epoch,E_LP_SC.data,E_JMAG_SC.data,LP_potentials[3],cross_coefficients_E)
    E_cal_SC = Struct(E_cal_SC,Epoch,None,'mV/m',['x - cal.','y - cal.','z - cal.'],'Calibrated E-field, from J-MAG')



# NOTE: This could use some work. Right now it does nothing and since it's made for a constant
# DC offset, which we might not have use of anymore, it may be unneccesary. 
"""     APPLY CALIBRATION/DC OFFSET
========================================================================
========================================================================
"""
"""
# Calculate calibrated differentials from E-field
# ================================================
LP_diffs_cal = np.matmul(rpwi_data.M_E2U,E_cal_SC.data)*1e-3
LP_diffs_cal = Struct(LP_diffs_cal,Epoch,Mask,'V',['P12','P23','P34'],'Calibrated LP differentials')


DC_offset = Struct(np.zeros((3,len(Epoch))),Epoch,Mask,'V',['Offset 1','Offset 2','Offset 3','Offset 4'],'DC offsets')


# Apply DC offsets to get calibrated potentials
# ================================================
LP_potentials_offset = Struct(np.zeros((4,len(Epoch))),Epoch,Mask,'V',['U1 + dU1', 'U2 + dU2', 'U3 + dU3', 'U4'],'Offset probe potentials')
LP_potentials_offset.data[3] = LP_potentials[3]
LP_potentials_offset.data[2] = LP_potentials[2] + DC_offset[2]
LP_potentials_offset.data[1] = LP_potentials[1] + DC_offset[1]
LP_potentials_offset.data[0] = LP_potentials[0] + DC_offset[0]


# Calculate new differentials using calibrated potentials
# ================================================
LP_diffs_offset = Struct(np.zeros((4,len(Epoch))),Epoch,Mask,'V',['P12','P23','P34','P4 (Single ended)'],'Offset probe differnetials')
LP_diffs_offset.data[3] = LP_diffs[3]
LP_diffs_offset.data[2] = LP_potentials_offset.data[2] - LP_potentials_offset.data[3]
LP_diffs_offset.data[1] = LP_potentials_offset.data[1] - LP_potentials_offset.data[2]
LP_diffs_offset.data[0] = LP_potentials_offset.data[0] - LP_potentials_offset.data[1]


# Calculate new E-field using offset differentials
# ================================================
E_LP_LPs_offset = Struct(np.zeros((3,len(Epoch))),Epoch,Mask,'mV/m',['P12', 'P23', 'P34'],'Offset E-field in LP differentials')
E_LP_LPs_offset.data[0] = 1e3*LP_diffs_offset.data[0]/np.linalg.norm(rpwi_data.LP12_distance)
E_LP_LPs_offset.data[1] = 1e3*LP_diffs_offset.data[1]/np.linalg.norm(rpwi_data.LP23_distance)
E_LP_LPs_offset.data[2] = 1e3*LP_diffs_offset.data[2]/np.linalg.norm(rpwi_data.LP34_distance)
"""



"""     PLOT (using matplotlib)
========================================================================
========================================================================
"""

# Toggle division of data into chunks.
# Plots new figures for each chunk, which might increse plotting performance.
# Alternative is to shorten specified time_period
# ================================================================
auto_chunking = False

# Limit size of chunks. I got 1024*512 (=524 288) from Ilona as a benchmark size, 
# but I found that 5e6 works well for me
chunk_size_limit = 5e6 
if auto_chunking:
    from math import ceil
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
        window_title = time_period[0].split('T')[0]
    else:
        window_title = time_period[0].split('T')[0] + ' Chunk ' + str(chunk_i+1) + ' (of ' + str(number_of_data_chunks) + ')'
    

    # Plot probe potentials
    # ================================================
    if False:
        fig, axes = plt.subplots(2,1,sharex=True)
        LP_potentials.plot(axes[0],None,chunk,4,plot_styles['lines'])
        LP_potentials_offset.plot(axes[1],None,chunk,4,plot_styles['lines'])
        fig.canvas.manager.set_window_title("LP potentials " + window_title)


    # Plot probe differentials [V]
    # ================================================
    if False:
        fig, axes = plt.subplots(3,1,sharex=True)
        LP_diffs.plot(axes[0],range(3),chunk,4,plot_styles['lines'])
        LP_diffs_offset.plot(axes[1],range(3),chunk,4,plot_styles['lines'])
        LP_diffs.plot(axes[2],range(3,4),chunk,4,plot_styles['lines'])
        fig.canvas.manager.set_window_title("LP differentials " + window_title)


    # Plot electric field strength in probe differentials
    # ================================================
    if False:
        fig, axes = plt.subplots(3,1,sharex=True)
        E_LP_LPs.plot(axes[0],None,chunk,4,plot_styles['lines'])
        E_LP_LPs_offset.plot(axes[1],None,chunk,4,plot_styles['lines'])
        LP_diffs.plot(axes[2],range(3,4),chunk,4,plot_styles['lines'])
        fig.canvas.manager.set_window_title("EFs " + window_title)


    # Plot electric field strenth in x,y,z-axes
    # ================================================
    if False:
        fig, axes = plt.subplots(3,1,sharex=True)
        E_LP_SC.plot(axes[0],None,chunk,4,plot_styles['lines'])
        E_cal_SC.plot(axes[1],None,chunk,4,plot_styles['lines'])
        LP_diffs.plot(axes[2],range(3,4),chunk,4,plot_styles['lines'])
        fig.canvas.manager.set_window_title("EFs in XYZ " + window_title)
    

    # IGRF model GSM
    # ========================================================================   
    if False:
        fig, axes = plt.subplots(2,1,sharex=True)
        B_IGRF_GSM.plot(axes[0])
        E_IGRF_GSM.plot(axes[1])

    
    # J-MAG B-field
    # ========================================================================  
    if True:
        fig, axes = plt.subplots(2,1,sharex=True)
        B_JMAG_SC.plot(axes[0])
        LP_potentials.plot(axes[1],range(3,4))



    # OMNI2 data
    # ========================================================================  
    if False:
        fig, axes = plt.subplots(2,1,sharex=True)
        B_OMNI_GSE.plot(axes[0])
        SW_OMNI.plot(axes[1],range(1))



    # J-MAG E-field
    # ========================================================================  
    if False:
        # E_LP_SC.data[0,:] = E_LP_SC.data[0,:] + 175
        # E_LP_SC.data[1,:] = E_LP_SC.data[1,:] + 30
        # E_LP_SC.data[2,:] = E_LP_SC.data[2,:] + 20
        fig, axes = plt.subplots(3,1,sharex=True)
        E_JMAG_SC.plot(axes[0],range(1,3))
        E_LP_SC.plot(axes[1],range(1,3))
        LP_potentials.plot(axes[2],range(3,4))



    # Fitted E-field
    # ========================================================================  
    if True:
        fig, axes = plt.subplots(3,1,sharex=True)
        E_LP_SC.plot(axes[0],range(1))
        E_cal_SC.plot(axes[0],range(1))
        axes[0].set_title('E-field fit to -v x B and V_sc compared to observed E, SC axis x')

        E_LP_SC.plot(axes[1],range(1,2))
        E_cal_SC.plot(axes[1],range(1,2))
        axes[1].set_title('E-field fit to -v x B and V_sc compared to observed E, SC axis y')

        E_LP_SC.plot(axes[2],range(2,3))
        E_cal_SC.plot(axes[2],range(2,3))
        axes[2].set_title('E-field fit to -v x B and V_sc compared to observed E, SC axis z')



    # Fitting coefficients
    # ========================================================================  
    if False:
        log_coefficients = True
        if log_coefficients:
            cal_coeff_a = abs(cal_coeff_a)
            cal_coeff_b = abs(cal_coeff_b)
            cal_coeff_c = abs(cal_coeff_c)


        for i in range(3):
            fig, axes = plt.subplots(3,1,sharex=True)
            E_LP_SC.plot(axes[0],range(i,i+1))
            E_cal_SC.plot(axes[0],range(i,i+1))
            axes[0].set_title('E-field fit to -v x B and V_sc compared to observed E, SC axis ' + ['x', 'y', 'z'][i])
            LP_potentials.plot(axes[2],range(3,4))
            if cross_coefficients_E:
                axes[1].plot(coeff_Epoc,cal_coeff_a.T[0][i], label = ['a_xx','a_yx','a_zx'][i], **plot_styles['lines'])
                axes[1].plot(coeff_Epoc,cal_coeff_a.T[1][i], label = ['a_xy','a_yy','a_zy'][i], **plot_styles['lines'])
                axes[1].plot(coeff_Epoc,cal_coeff_a.T[2][i], label = ['a_xz','a_yz','a_zz'][i], **plot_styles['lines'])
            else:
                axes[1].plot(coeff_Epoc,cal_coeff_a.T[i][i], label = ['a_xx','a_yy','a_zz'][i], **plot_styles['lines'])
            axes[1].plot(coeff_Epoc,cal_coeff_b.T[i],label = 'b_x', **plot_styles['lines'])
            axes[1].plot(coeff_Epoc,cal_coeff_c.T[i],label = 'c_x', **plot_styles['lines'])
            axes[1].legend(loc='upper right')
            axes[1].grid(True)
              
            if log_coefficients:
                axes[1].set_yscale('log')
        

    # Compare with cyclotron frequency
    # ======================================================================== 
    if False:
        time_period_cyclotron = ['2024-08-20T21:07:59','2024-08-20T23:45']
        included_data_cyclotron = (Epoch_LEGA > spice.unitim(spice.utc2et(time_period_cyclotron[0]), "ET", "TT")*1e9) & (Epoch_LEGA < spice.unitim(spice.utc2et(time_period_cyclotron[1]), "ET", "TT")*1e9)
        e = 1.6e-19
        m_e = 9.109e-31
        freq_ce = Struct(np.array([(e/(2*np.pi*m_e))*np.linalg.norm(B_IGRF_SC.data[:,included_data_cyclotron],axis=0)*1e-9]),Epoch_LEGA[included_data_cyclotron],None,'Hz',['f_ce'],'Electron gyrofrequency from IGRF B-field')
        
        # Save
        from scipy.io import savemat
        savemat(rootDir + "/fsweep/freq_ce.mat",{'freq_ce': freq_ce.data,'time':freq_ce.time})

        # Plot
        plt.figure()
        freq_ce.plot(plt.gca())
        plt.yscale('log')
        plt.grid(True, which='both')


    # This chunk's plots
    # ================================================
    plt.show()

# print(np.mean(cal_coeff_a))