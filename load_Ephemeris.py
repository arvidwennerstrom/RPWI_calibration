import scipy, numpy as np
from E_calibr_support import datenum_to_tt2000

def load_Ephemeris(rootDir,Epoch):

     # Get Epoch for ephemeris during LEGA
    Epoch_ephm = []
    for t in scipy.io.loadmat(rootDir + '/Ephemeris/datenum_LEGA.mat')['datenum'][0]:
        Epoch_ephm.append(datenum_to_tt2000(t))
    Epoch_ephm = 1e9*np.array(Epoch_ephm)
    

    # Truncate data to only use specified time period
    included_data_LEGA = (Epoch_ephm > Epoch[0]) & (Epoch_ephm < Epoch[-1])
    Epoch_ephm = Epoch_ephm[included_data_LEGA]

  
    # Load ephemeris data
    LEGA_ephemeris = scipy.io.loadmat(rootDir + '/Ephemeris/LEGA_ephemeris 2024-11-12 10_30_16.mat')
    EARTH_SC = LEGA_ephemeris['eph']['EARTH_SC'][0][0][:,included_data_LEGA]
    pos_ephem_GSM = LEGA_ephemeris['eph']['juice_GSM'][0][0][:,included_data_LEGA]
    v_ephem_GSM = LEGA_ephemeris['eph']['juice_V_GSM'][0][0][:,included_data_LEGA] 


    return EARTH_SC, pos_ephem_GSM, v_ephem_GSM, Epoch_ephm