import numpy as np, spiceypy as spice
from E_calibr_support import Struct, year_doy_to_tt2000
from load_Ephemeris import load_Ephemeris


def load_omni_data(rootDir,Epoch):


    # Read OMNI2 data
    # ========================================================================  
    explanation = open(rootDir + '/OMNI2/omni2_explanation.txt',"r").read()
    data =  open(rootDir + '/OMNI2/omni2_data.txt',"r").read().split('\n')[:-1]

    t = []; B_GSM = []; SW_abs = []
    for idx,line in enumerate(data):
        clean_line = [val for val in line.split(' ') if val != '']
        t.append(year_doy_to_tt2000(int(clean_line[0]), int(clean_line[1]), int(clean_line[2])))
        B_GSM.append(clean_line[3:6])
        SW_abs.append(float(clean_line[6]))


    # Convert data to np.arrays, for easier handling
    t = np.array(t)
    B_GSM = np.array(B_GSM).astype(float).T
    SW_abs = np.array(SW_abs)
    


    # Temporally shift Epoch to account for difference in Juice's and L1's 
    # radial distance from the Sun.
    # ======================================================================== 
    EARTH_SC, juice_pos_GSM_RE, v_ephem_GSM, Epoch_ephm = load_Ephemeris(rootDir,Epoch)
    
    # Get Juice's position in km in GSM
    R_E = 6371.008366666666 # [km]
    juice_pos_GSM_km = juice_pos_GSM_RE*R_E
    
    # Get the Sun's position in km in GSM
    AU = 1.49597871*1e8 # [km]
    Sun_pos_GSM_km = np.zeros_like(juice_pos_GSM_km)
    Sun_pos_GSM_km[0] = Sun_pos_GSM_km[0] + AU

    Sun_2_juice = juice_pos_GSM_km - Sun_pos_GSM_km
    r_juice = np.average(np.linalg.norm(Sun_2_juice, axis = 0)) # Radial distance of Juice from the Sun [km]


    r_L1 = 1.481*1e8 # Radial distance of L1 from the Sun [km]
    
    r_diff = r_juice-r_L1 # Radial difference between L1 and Juice

    

    # NOTE: Velocity is set to 340 km/s, because it approximates the SW of the period 04:00-05:30 on August 23.
    # Maybe improve this
    v = 340 # [km/s]


    # Shift the time of the omni data
    t_shift = r_diff/v # [seconds]
    Epoch = Epoch-t_shift*1e9


    # Only include data within specified time period
    # ======================================================================== 
    included_data = (t > Epoch[0]) & (t < Epoch[-1])
    t = t[included_data]
    B_GSM = B_GSM[:,included_data]
    SW_abs = SW_abs[included_data]


    # Remove data where value saturates
    # ========================================================================  
    B_GSM[B_GSM == 999.9] = np.nan
    SW_abs[SW_abs == 9999] = np.nan



    # Create 3D SW-data. SW_abs is just magnitude, from Sun to Earth, along -x axis in GSE
    # ========================================================================  
    SW_GSM = np.zeros_like(B_GSM)
    SW_GSM[0] = -SW_abs




    # Create Struct-type objects and return 
    # ========================================================================  
    B = Struct(B_GSM,t,None,'nT',['X_GSE','Y_GSE','Z_GSE'],'ACE magnetic field in GSE')
    SW = Struct(SW_GSM,t,None,'km/s',['Sw_GSEx', 'V_sw_y', 'V_sw_z'],'SW plasma speed in GSE')
    return B, SW



# # Extrapolate data to fit the higher resoultion measured data
# # ========================================================================  
# B_high_res = np.zeros((3,len(Epoch)))
# SW_high_res = np.zeros((1,len(Epoch)))


# for idx,b in enumerate(B_GSM):
#     B_high_res[idx] = np.interp(Epoch,t,B_GSM[idx])

# SW_high_res[0,:] = np.interp(Epoch,t,SW_abs)