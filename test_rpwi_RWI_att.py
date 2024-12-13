import numpy as np
import spiceypy as spice
import datetime
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

# spice.furnsh('/Users/morooka/DATA/SPICE/mk_files/mk_juice/juice_ops.tm')

# Input overall directory path, where folders "datasets", 
# "data_created" and "spice" exist
# ========================================================================
rootDir = "C:/Users/1/Onedrive - KTH/MEX/IRF" # My desktop
rootDir = "C:/Users/arvwe/Onedrive - KTH/MEX/IRF" # My laptop
# rootDir = "C:/Users/arvidwen/Onedrive - KTH/MEX/IRF" # KTH computers

spice.furnsh((rootDir + '/spice/JUICE/kernels/juice_ops.tm'))


tseries = pd.date_range(start='2024-08-19 00:00:00',end='2024-08-28 00:00:00',freq='60s')
et      = spice.str2et(tseries.strftime("%Y-%m-%d %H:%M:%S"))


# --- JUICE POSITION in GSE ---
JUICE_GSE, ltime = spice.spkpos('JUICE',et,'GSE','LT+S','EARTH')

# SUN/Earth pointing in spacecraft coordinates
SUN_POS,  ltime = spice.spkpos('SUN',et,'JUICE_SPACECRAFT','LT+S','JUICE');
norm = np.array([np.linalg.norm(SUN_POS,axis=1)]).T * np.array([1,1,1])
SUN_dir = np.divide(SUN_POS,norm)

Earth_POS,  ltime = spice.spkpos('EARTH',et,'JUICE_SPACECRAFT','LT+S','JUICE');   
norm = np.array([np.linalg.norm(Earth_POS,axis=1)]).T * np.array([1,1,1])
Earth_dir = np.divide(Earth_POS,norm)

fig, ax = plt.subplots(2)
for ii in np.arange(0,3):
	ax[1].plot(tseries, SUN_dir[:,ii])
	ax[0].plot(tseries, Earth_dir[:,ii])


# SUN/Earth pointing in RWI coordinates
SUN_POS,  ltime = spice.spkpos('SUN',et,'JUICE_MAG_BOOM','LT+S','JUICE')
norm = np.array([np.linalg.norm(SUN_POS,axis=1)]).T * np.array([1,1,1])
SUN_dir = np.divide(SUN_POS,norm)

Earth_POS,  ltime = spice.spkpos('EARTH',et,'JUICE_MAG_BOOM','LT+S','JUICE');   
norm = np.array([np.linalg.norm(Earth_POS,axis=1)]).T * np.array([1,1,1])
Earth_dir = np.divide(Earth_POS,norm)

ll = np.arange(0,len(et))
matrix=np.zeros((3,3,len(et)))
for ii in ll:
	matrix[:,:,ii] = spice.pxform('JUICE_MAG_BOOM','JUICE_SPACECRAFT',et[ii])

# ========================================================================
print(matrix.T[0])
# ========================================================================

earthdir_new =np.zeros((len(et),3))
sundir_new   =np.zeros((len(et),3))
for ii in ll:
	earthdir_new[ii,:] = np.dot(matrix[:,:,ii], Earth_dir[ii,:])
	sundir_new[ii,:]   = np.dot(matrix[:,:,ii], SUN_dir[ii,:])

fig, ax = plt.subplots(2)
for ii in np.arange(0,3):
	ax[1].plot(tseries, earthdir_new[:,ii])
	ax[0].plot(tseries, sundir_new[:,ii])

spice.kclear()

plt.show()


