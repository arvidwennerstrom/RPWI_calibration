import numpy as np, spiceypy as spice
from E_calibr_support import Struct


def load_jmag_data(rootDir,Epoch):
    

    # Read J-MAG data
    # ========================================================================  
    data = open(rootDir + '/J-MAG/2024-08-23_fob_16Hz.dat',"r").read().split('\n')[:-1]

    
    t = []; B = np.zeros((4,len(data)))
    for idx,line in enumerate(data):
        vals = line.split('\t')
        

        # Read and convert time to tt2000
        tt2000 = spice.unitim(spice.str2et(vals[0]), 'ET', 'TT')*1e9
        t.append(tt2000)
        

        # Read B-field data
        B[:,idx] = vals[1:]

    t = np.array(t)


    # Convert from J-MAG science coordinates to SC coords. R-matrix given by J-MAG   
    R = np.array([
        [6.29320391*1e-1, -3.88573117*1e-1,  6.73028067*1e-1],
        [7.70695203*1e-17,  8.66025303*1e-1,  5.00000175*1e-1],
        [-7.77145961*1e-1, -3.14660306*1e-1,  5.45007382*1e-1]])
    B[0:3,:] = np.matmul(R,B[0:3,:])


    # Only include data within specified time period
    # ========================================================================  
    included_data = (t > Epoch[0]) & (t < Epoch[-1])
    t = t[included_data]
    B = B[:,included_data]



    # Extrapolate data to fit the higher resoultion measured data
    # ========================================================================  
    B_high_res = np.zeros((4,len(Epoch)))

    # Loop over x,y,z and extrapolate to match the size of Epoch
    for idx,b in enumerate(B):
        B_high_res[idx] = np.interp(Epoch,t,B[idx])



    # Create a Struct-type object and return 
    # ========================================================================  
    B_struct = Struct(B_high_res,Epoch,None,'nT',['Bx','By','Bz','|B|'],'J-MAG magnetic field in SC coordinates')
    return B_struct

