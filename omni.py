import numpy as np, spiceypy as spice
from E_calibr_support import Struct, year_doy_to_tt2000
from load_Ephemeris import load_Ephemeris


def load_omni_data(rootDir,Epoch):


    # Read OMNI2 data
    # ========================================================================  
    explanation = open(rootDir + '/OMNI2/omni2_explanation.txt',"r").read()
    data =  open(rootDir + '/OMNI2/omni2_data.txt',"r").read().split('\n')[:-1]

    t_orig = []; B_orig = []; SW_orig = []
    for idx,line in enumerate(data):
        clean_line = [val for val in line.split(' ') if val != '']
        t_orig.append(year_doy_to_tt2000(int(clean_line[0]), int(clean_line[1]), int(clean_line[2])))
        B_orig.append(clean_line[3:6])
        SW_orig.append(float(clean_line[6]))


    # Convert data to np.arrays, for easier handling
    t_orig = np.array(t_orig)
    B_orig = np.array(B_orig).astype(float).T
    SW_orig = np.array(SW_orig)
    


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
    included_data = (t_orig > Epoch[0]) & (t_orig < Epoch[-1])
    t_included = t_orig[included_data]
    B_included = B_orig[:,included_data]
    SW_included = SW_orig[included_data]


    
    # Extrapolate data to fit the higher resoultion measured data
    # ========================================================================  
    B_high_res = np.zeros((3,len(Epoch)))
    SW_high_res = np.zeros((1,len(Epoch)))

    # outside_range = (Epoch < t_included[0]) | (Epoch > t_included[-1])


    for idx,b in enumerate(B_included):
        B_high_res[idx] = np.interp(Epoch,t_included,B_included[idx])
        # B_high_res[idx][outside_range] = np.nan

    SW_high_res[0,:] = np.interp(Epoch,t_included,SW_included)
    # SW_high_res[outside_range] = np.nan


    # Remove data where value saturates
    # ========================================================================  
    B_high_res[B_high_res == 999.9] = np.nan
    SW_high_res[SW_high_res == 9999] = np.nan



    # Create Struct-type objects and return 
    # ========================================================================  
    B = Struct(B_high_res,Epoch,None,'nT',['X_GSE','Y_GSE','Z_GSE'],'ACE magnetic field in GSE')
    SW = Struct(SW_high_res,Epoch,None,'km/s',['V_sw_x', 'V_sw_y', 'V_sw_z'],'SW plasma speed in GSE')
    return B, SW


