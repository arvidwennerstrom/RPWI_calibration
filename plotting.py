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
import numpy as np, scipy, matplotlib.pyplot as plt
from datetime import datetime, timedelta
from rpwi_data import rpwi_data, coeffs_TM2voltage  
from plotting_support import Struct, datenum_to_tt2000, rotation_matrix



"""     INPUT OPTIONS AND SETTINGS FOR PLOTTING
========================================================================           
========================================================================
"""
# Input overall directory path, where folders "datasets", 
# "data_created" and "spice" exist
# ========================================================================
rootDir = "C:/Users/1/Onedrive - KTH/MEX/IRF" # My desktop
rootDir = "C:/Users/arvwe/Onedrive - KTH/MEX/IRF" # My laptop
# rootDir = "C:/Users/arvidwen/Onedrive - KTH/MEX/IRF" # KTH computers


# Input date and time (optional) of the data to plot.
# ========================================================================
# Specify on format: "YYYY-MM-DDTHH:MM:SS.ns*". 
# Date is mandatory, adding T is optional: "2024-08-20" and 2024-08-20T10" are both valid. 
# Time precision is optional: "*T08" and "*T08.30:20.500" are both valid. 
time_period = ['2024-08-20T20', '2024-08-21T01']

start_date = datetime.strptime(time_period[0].split('T')[0], "%Y-%m-%d")
end_date = datetime.strptime(time_period[1].split('T')[0], "%Y-%m-%d")
date_list = [(start_date + timedelta(days=x)).strftime("%Y%m%d") for x in range((end_date - start_date).days + 1)]


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
LP_potentials = Struct(np.zeros((4,len(Epoch))),Epoch,Mask,'V',['P01 (P12+P23+P34+P04)', 'P02 (P23+P34+P04)', 'P03 (P34+P04)', 'P04'],'Probe potentials SID1')
LP_potentials.data[3] = LP_diffs[3]
LP_potentials.data[2] = LP_potentials.data[3] + LP_diffs.data[2]
LP_potentials.data[1] = LP_potentials.data[2] + LP_diffs.data[1]
LP_potentials.data[0] = LP_potentials.data[1] + LP_diffs.data[0]



"""     ELECTRIC FIELD
======================================================================== 
========================================================================
"""
# M is transformation matrix in: U = M*E, where U is voltage and E is electric field
M = [rpwi_data.LP12_distance,
    rpwi_data.LP23_distance,
    rpwi_data.LP34_distance]
M_inv = np.linalg.inv(M)


# Get E-field between probes
# ================================================ 
# 1e3's are used to give EF in mV/m, rather than V/m
E_LP_LPs = Struct(np.zeros((3,len(Epoch))),Epoch,Mask,'mV/m',['P12', 'P23', 'P34'],'E-field in LP differentials')
for i in range(3):
    E_LP_LPs.data[i] = 1e3*LP_diffs[i]/np.linalg.norm(M[i]) 
    

# Calculate E-field in SC x-,y- and z-axes
# ================================================
E_LP_SC = Struct(1e3*np.matmul(M_inv,LP_diffs.data[:3]),Epoch,Mask,'mV/m',['x','y','z'],'LP E-field in SC coords.')



"""     CALCULATE DC OFFSET 
========================================================================
========================================================================
"""
""" UNUSED
# Sweep data
# NOTE: This should use regular expr.
# ========================================================================
# if date == '240820':
#     mat_data = scipy.io.loadmat(rootDir + '/sweep_data/lpstruc_20' + date + '_2_20241002.mat')
#     t_sweep = mat_data['lpstruc']['t_sweep'][0][0]
#     Usc = mat_data['lpstruc']['Usc'][0][0]
#     Ni = abs(mat_data['lpstruc']['Ni'][0][0])

#     Epoch_sweep = [datenum_to_tt2000(t[0])*1e9 for t in t_sweep]  
"""


# IGRF B-field (using geopack)
# ========================================================================
# Documentation on geopack: (https://github.com/tsssss/geopack/?tab=readme-ov-file)
from geopack import geopack

# Closest approach @ 2024-08-20T21:56:14
Unix_CA = 1724190974 

# Set up geopack for correct time. Geomagnetic field is considered stationary for simplicty
geopack.recalc(Unix_CA) 

# Get Epoch for ephemeris during LEGA
Epoch_LEGA = []
for t in scipy.io.loadmat(rootDir + '/Ephemeris/datenum_LEGA.mat')['datenum'][0]:
    Epoch_LEGA.append(datenum_to_tt2000(t))
Epoch_LEGA = 1e9*np.array(Epoch_LEGA)


# Truncate data to only use specified time period
included_data_LEGA = (Epoch_LEGA > spice.unitim(spice.utc2et(time_period[0]), "ET", "TT")*1e9) & (Epoch_LEGA < spice.unitim(spice.utc2et(time_period[1]), "ET", "TT")*1e9)
Epoch_LEGA = Epoch_LEGA[included_data_LEGA]


# Load ephemeris data
LEGA_ephemeris = scipy.io.loadmat(rootDir + '/Ephemeris/LEGA_ephemeris 2024-11-12 10_30_16.mat')
EARTH_SC = Struct(LEGA_ephemeris['eph']['EARTH_SC'][0][0][:,included_data_LEGA],Epoch_LEGA,None,'',['x','y','z'],'Earth direction in SC coords.')
pos_ephem_GSE = Struct(LEGA_ephemeris['eph']['juice_GSE'][0][0][:,included_data_LEGA],Epoch,None,'R_E',['X','Y','Z'],"JUICE's position in GSE coordinates")
pos_ephem_GSM = Struct(LEGA_ephemeris['eph']['juice_GSM'][0][0][:,included_data_LEGA],Epoch,None,'R_E',['X','Y','Z'],"JUICE's position in GSM coordinates")
v_ephem_GSE = Struct(LEGA_ephemeris['eph']['juice_V_GSE'][0][0][:,included_data_LEGA],Epoch,None,'km/s',['X','Y','Z'],"JUICE's velocity in GSE coordinates")
v_ephem_GSM = Struct(LEGA_ephemeris['eph']['juice_V_GSM'][0][0][:,included_data_LEGA],Epoch,None,'km/s',['X','Y','Z'],"JUICE's velocity in GSM coordinates") 


# Naming convention: *typeofdata*_*obtainedfrom*_*coordinatesystem* 
# Example: B_IGRF_GSM is magnetic field data, from IGRF model, in GSM coordinates
B_IGRF_GSM = Struct(np.zeros((4,len(Epoch_LEGA))),Epoch_LEGA,None,'nT',['X (from Earth to Sun)','Y (in magnetic equatorial plane)',"Z (aligned with Earth's mag. dipole axis)",'Magnitude'],'IGRF B-field in GSM')
B_IGRF_SC = Struct(np.zeros((4,len(Epoch_LEGA))),Epoch_LEGA,None,'nT',['x (away from Sun)','y','z','Magnitude'],'IGRF B-field in SC coords.')
E_IGRF_GSM = Struct(np.zeros((4,len(Epoch_LEGA))),Epoch_LEGA,None,'mV/m',['X','Y','Z','Magnitude'],'IGRF E-field in GSM')
E_IGRF_SC = Struct(np.zeros((4,len(Epoch_LEGA))),Epoch_LEGA,None,'mV/m',['x (away from Sun)','y','z','Magnitude'],'IGRF E-field in SC coords.')
juice_V_SC = Struct(np.zeros((3,len(Epoch_LEGA))),Epoch_LEGA,None,'km/s',['x (away from Sun)','y','z'],'Juice velocity in SC coords.')


# List of rotation matrices.
R_GSM2SC = []


# Make calculations separately for every step of ephemeris data
for i in range(len(Epoch_LEGA)):
    
    # Unit vectors for direction of Earth in Juice coordinates and direction
    # of Juice in GSM coordinates. These are pointing along the same axis,
    # but in opposite directions, which is why -pos_ephem_GSM is needed.  
    u_SC = EARTH_SC[:,i]/np.linalg.norm(EARTH_SC[:,i])
    u_GSM = -pos_ephem_GSM[:,i]/np.linalg.norm(pos_ephem_GSM[:,i])

    # Rotation matrix from GSM to SC coords.
    R_GSM2SC.append(rotation_matrix(u_GSM,u_SC))

    # B-field 
    # Get B-field from IGRF.
    # Perform coordinate transformations.
    xgsm,ygsm,zgsm = pos_ephem_GSM[:,i]
    bgsm = np.array(geopack.igrf_gsm(xgsm,ygsm,zgsm))
    B_IGRF_GSM.data[0:3,i] = bgsm
    B_IGRF_GSM.data[3,i] = np.linalg.norm(bgsm)
    
    bsc = np.dot(R_GSM2SC[i],bgsm)
    B_IGRF_SC.data[0:3,i] = bsc
    B_IGRF_SC.data[3,i] = np.linalg.norm(bsc)

    # Velocity
    # Get V from ephemeris.
    # Perform coordinate transformations.
    vgsm = v_ephem_GSM[:,i]
    vsc = np.dot(R_GSM2SC[i],vgsm)
    juice_V_SC.data[:,i] = vsc


    # E-field
    # Get E-field from: E = -v x B
    # v is in [km/s], B in [nT] and desired output is E in [mV/m]. Scaling with 1e-3 ensures this.
    Esc = np.cross(-vsc,bsc)*1e-3
    E_IGRF_SC.data[0:3,i] = Esc
    E_IGRF_SC.data[3,i] = np.linalg.norm(Esc)

    Egsm = np.cross(-vgsm,bgsm)*1e-3
    E_IGRF_GSM.data[0:3,i] = Egsm
    E_IGRF_GSM.data[3,i] = np.linalg.norm(Egsm)
    


# Calculate DC offset
# ================================================
# DC offsets is a np.array corresponding to data values, 
# to allow changes in offset with time. 
# NOTE: P4 has an offset to in case it is needed in the future, though it is zero atm 
DC_offset = Struct(np.zeros((4,len(Epoch))),Epoch,Mask,'V',['Offset 1','Offset 2','Offset 3','Offset 4'],'DC offsets')

"""
Anders' idea: 
    E_LPs = A*E_vxB + B*P4 + C
    E_LPs is measured E field by LP's, E_vxB is "true" E-field from E = -v x B,
    U_SC is spacecraft potential. 
    A, B and C are constants to be found, by Least square or equal
Rewrite as optimization problem:
    minimize: sum(E_LPs - (A*E_vxB + B*P4 + C))^2


Possibly add D*(1/P4). I don't know if it makes sense, but I think it will work well 
"""

# Create coefficients. They are of length 4 to correspond to each LP channel: P12, P23, P34 and P4
A = np.zeros(4)
B = np.zeros(4)
C = np.zeros(4)

# # For testing 1/P4
# D = np.zeros(4)


# Extend vxB-data
E_IGRF_SC_long = Struct(np.zeros((3,len(Epoch))),Epoch,None,'mv/m',['x','y','z'],'E-field in SC coords.')
for i in range(3):
    E_IGRF_SC_long.data[i] = np.interp(Epoch,Epoch_LEGA,E_IGRF_SC[i])

# P4_inv = LP_potentials[3]**(-1)


E_LP_SC_offset = Struct(np.zeros((3,len(Epoch))),Epoch,Mask,'mV/m',['x - offset','y - offset','z - offset'],'Offset LP E-field in SC coords')
# Offset_Test = np.zeros((3,len(Epoch)))

# Find coefficient values using least square fitting
for i in range(3):
    # Create the design matrix X with shape (len(Epoch_LEGA), 3)
    # The columns represent [E_vxB, P4, 1] for the coefficients A, B, C respectively
    X = np.vstack([E_IGRF_SC_long[i],LP_potentials[3],np.ones(len(Epoch))]).T

    # Solve using least square sum
    coeffs, residual, rank, s = np.linalg.lstsq(X,E_LP_SC[i],rcond=None)
    A[i],B[i],C[i] = coeffs

    E_LP_SC_offset.data[i] = A[i]*E_IGRF_SC_long[i] + B[i]*LP_potentials[3] + C[i] 



fig, axes = plt.subplots(3,1,sharex=True)
E_LP_SC.plot(axes[0],range(1))
E_LP_SC_offset.plot(axes[0],range(1))
axes[0].set_title('Comparison between observed and offset E-field, x')
axes[0].legend(loc='upper right')

E_LP_SC.plot(axes[1],range(1,2))
E_LP_SC_offset.plot(axes[1],range(1,2))
axes[1].set_title('Comparison between observed and offset E-field, y')
axes[1].legend(loc='upper right')

E_LP_SC.plot(axes[2],range(2,3))
E_LP_SC_offset.plot(axes[2],range(2,3))
axes[2].set_title('Comparison between observed and offset E-field, z')
axes[2].legend(loc='upper right')



"""     APPLY DC OFFSET 
========================================================================
========================================================================
"""

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

# E_LP_SC_offset = Struct(1e3*np.matmul(M_inv,LP_diffs_offset.data[:3]),Epoch,Mask,'mV/m',['x','y','z'],'Offset LP E-field in SC coords.')


"""     DOWNSAMPLING OF OBS-DATA
======================================================================== 
========================================================================
"""
# Downsample observed data and convert it to GSM
E_obs_SC_short = np.zeros((3,len(Epoch_LEGA)))
for axis in range(3):
    E_obs_SC_short[axis] = np.interp(Epoch_LEGA,Epoch,E_LP_SC.data[axis])

E_obs_GSM_short = np.zeros((3,len(Epoch_LEGA)))
for i in range(len(Epoch_LEGA)):
    E_obs_GSM_short[:,i] = np.dot(np.linalg.inv(R_GSM2SC[i]),E_obs_SC_short[:,i])



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
    

    # Plot SC potential at each probe
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
        E_LP_SC_offset.plot(axes[1],None,chunk,4,plot_styles['lines'])
        LP_diffs.plot(axes[2],range(3,4),chunk,4,plot_styles['lines'])
        fig.canvas.manager.set_window_title("EFs in XYZ " + window_title)
    

    # IGRF model GSM
    # ========================================================================   
    if False:
        fig, axes = plt.subplots(4,1,sharex=True)
        B_IGRF_GSM.plot(axes[0])
        E_IGRF_GSM.plot(axes[1])

        axes[2].plot(Epoch_LEGA,E_obs_GSM_short[0],label='X')
        axes[2].plot(Epoch_LEGA,E_obs_GSM_short[1],label='Y')
        axes[2].plot(Epoch_LEGA,E_obs_GSM_short[2],label='Z')
        axes[2].set_ylabel('mv/m')
        axes[2].legend(loc='upper right')
        axes[2].grid(True)
        axes[2].set_title('JUICE observed E-field in GSM')

        # Position
        axes[3].plot(Epoch_LEGA,pos_ephem_GSM[0],label='X')
        axes[3].plot(Epoch_LEGA,pos_ephem_GSM[1],label='Y')
        axes[3].plot(Epoch_LEGA,pos_ephem_GSM[2],label='Z')
        axes[3].set_ylabel('R_E')
        axes[3].legend(loc='upper right')
        axes[3].grid(True)
        axes[3].set_title('JUICE position in GSM')


    # Downsampling visualization
    # ========================================================================  
    if False:
        fig, axes = plt.subplots(2,1,sharex=True)
        E_LP_SC.plot(axes[0])
        axes[1].plot(Epoch_LEGA,E_obs_SC_short[0],label='X')
        axes[1].plot(Epoch_LEGA,E_obs_SC_short[1],label='Y')
        axes[1].plot(Epoch_LEGA,E_obs_SC_short[2],label='Z')
        axes[1].legend(loc='upper right')
        axes[1].set_ylabel('mv/m')
        axes[1].grid(True)
        axes[1].set_title('Downsampled E-field')
    

    # Compare with cyclotron frequency
    # ======================================================================== 
    if False:
        time_period_cyclotron = ['2024-08-20T20:40','2024-08-20T23']
        included_data_cyclotron = (Epoch_LEGA > spice.unitim(spice.utc2et(time_period_cyclotron[0]), "ET", "TT")*1e9) & (Epoch_LEGA < spice.unitim(spice.utc2et(time_period_cyclotron[1]), "ET", "TT")*1e9)
        e = 1.6e-19
        m_e = 9.109e-31
        freq_ce = Struct(np.array([(e/(2*np.pi*m_e))*np.linalg.norm(B_IGRF_SC.data[:,included_data_cyclotron],axis=0)*1e-9]),Epoch_LEGA[included_data_cyclotron],None,'Hz',['f_ce'],'Electron gyrofrequency from IGRF B-field')
        plt.figure()
        freq_ce.plot(plt.gca())
        plt.yscale('log')
        plt.grid(True, which='both')




    # fig, axes = plt.subplots(4,1,sharex=True)
    # B_IGRF_SC.plot(axes[0])
    # E_IGRF_SC.plot(axes[1])
    # juice_V_SC.plot(axes[2])
    # axes[3].plot(Epoch_LEGA,EARTH_SC[0,:],label='x')
    # axes[3].plot(Epoch_LEGA,EARTH_SC[1,:],label='y')
    # axes[3].plot(Epoch_LEGA,EARTH_SC[2,:],label='z')
    # axes[3].legend(loc='upper right')
    # axes[3].grid(True)
    # axes[3].set_title('Earth direction in SC coords.')


    # IGRF model B-field
    # ========================================================================        
    # fig, axes = plt.subplots(3,1,sharex=True)
    # B_IGRF_GSM.plot(axes[0])
    # B_IGRF_SC.plot(axes[1])

    # axes[2].plot(Epoch_LEGA,EARTH_SC[0,:],label='X')
    # axes[2].plot(Epoch_LEGA,EARTH_SC[1,:],label='Y')
    # axes[2].plot(Epoch_LEGA,EARTH_SC[2,:],label='Z')
    # axes[2].legend(loc='upper right')
    # axes[2].grid(True)
    # axes[2].set_title('Earth direction in JUICE axes')


    # # Comparison between measured and IGRF E-field
    # fig, axes = plt.subplots(2,1,sharex=True)
    # E_LP_SC.plot(axes[0],chunk)
    # E_IGRF_SC.plot(axes[1])


    # fig, axes = plt.subplots(2,1,sharex=True)
    # B_IGRF_GSM.plot(axes[0])
    # E_IGRF_GSM.plot(axes[1])



    # # Verification purposes
    # fig, axes = plt.subplots(2,1,sharex=True)
    # axes[0].plot(pos_ephem_GSE[0],label='X')
    # axes[0].plot(pos_ephem_GSE[1],label='Y')
    # axes[0].plot(pos_ephem_GSE[2],label='Z')
    # axes[0].legend(loc='upper right')
    # axes[0].grid(True)
    # axes[0].set_title('Juice velocity in GSE')



    # This chunk's plots
    # ================================================
    plt.show()

