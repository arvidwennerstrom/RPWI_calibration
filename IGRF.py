
def IGRF_B_field(rootDir,Epoch):
    import numpy as np
    from geopack import geopack
    from E_calibr_support import Struct, rotation_matrix, tt2000_to_unix
    from load_Ephemeris import load_Ephemeris
    # IGRF B-field using geopack
    # Documentation on geopack: (https://github.com/tsssss/geopack/?tab=readme-ov-file)


    EARTH_SC, pos_ephem_GSM, v_ephem_GSM, Epoch_ephm = load_Ephemeris(rootDir,Epoch)


    # List of rotation matrices.
    R_GSM2SC = np.zeros((len(Epoch_ephm),3,3))


    # Naming convention: *typeofdata*_*obtainedfrom*_*coordinatesystem* 
    # Example: B_GSM is magnetic field data, from IGRF model, in GSM coordinates
    B_GSM = np.zeros((len(Epoch_ephm),3))
    v_GSM = v_ephem_GSM.T


    # Make calculations separately for every step of ephemeris data
    for i in range(len(Epoch_ephm)):
        
        # Recalc IGRF for each time step
        unix = tt2000_to_unix(Epoch_ephm[i])
        geopack.recalc(unix)


        # Get B-field from IGRF and calculate E-field, in GSM
        # ================================================
        # B-field from IGRF
        xgsm,ygsm,zgsm = pos_ephem_GSM[:,i]
        bgsm = np.array(geopack.igrf_gsm(xgsm,ygsm,zgsm))
        B_GSM[i,:] = bgsm


        # Calculate transformation matrix from GSM to SC coords.
        # ================================================
        # Unit vectors for direction of Earth in Juice coordinates and direction
        # of Juice in GSM coordinates. These are pointing along the same axis,
        # but in opposite directions, which is why -pos_ephem_GSM is needed. 
        # NOTE: Am I just confused, or is this wrong? So "pos_ephem_GS" and not "-pos_ephem_GSM"?
        u_SC = EARTH_SC[:,i]/np.linalg.norm(EARTH_SC[:,i])
        u_GSM = pos_ephem_GSM[:,i]/np.linalg.norm(pos_ephem_GSM[:,i])
        
        # Rotation matrix from GSM to SC coords.
        R = rotation_matrix(u_GSM,u_SC)
        R_GSM2SC[i,:,:] = R


    # Get E-field in GSM
    # ================================================
    E_GSM = np.cross(-v_GSM,B_GSM)


    # Convert to SC coordinates
    # ================================================
    B_SC = np.matmul(R_GSM2SC,B_GSM[...,None])[..., 0]
    v_SC = np.matmul(R_GSM2SC,v_GSM[...,None])[..., 0]
    E_SC = np.matmul(R_GSM2SC,E_GSM[...,None])[..., 0]



    # Extend length of data to match resolution of LP measurements
    # ================================================
    B_GSM_long = np.zeros((len(Epoch),4))
    B_SC_long = np.zeros((len(Epoch),4))
    E_GSM_long = np.zeros((len(Epoch),4))
    E_SC_long = np.zeros((len(Epoch),4))
    v_SC_long = np.zeros((len(Epoch),4))


    for i in range(3):
        B_GSM_long[:,i] = np.interp(Epoch,Epoch_ephm,B_GSM[:,i])
        B_SC_long[:,i] = np.interp(Epoch,Epoch_ephm,B_SC[:,i])
        E_GSM_long[:,i] = np.interp(Epoch,Epoch_ephm,E_GSM[:,i])
        E_SC_long[:,i] = np.interp(Epoch,Epoch_ephm,E_SC[:,i])
        v_SC_long[:,i] = np.interp(Epoch,Epoch_ephm,v_SC[:,i])



    # Add the norm of x,y and z components
    # ================================================
    B_GSM_long[:,3] = np.linalg.norm(B_GSM_long[:,0:3],axis=1)
    B_SC_long[:,3] = np.linalg.norm(B_SC_long[:,0:3],axis=1)
    E_GSM_long[:,3] = np.linalg.norm(E_GSM_long[:,0:3],axis=1)
    E_SC_long[:,3] = np.linalg.norm(E_SC_long[:,0:3],axis=1)
    v_SC_long[:,3] = np.linalg.norm(v_SC_long[:,0:3],axis=1)



    # Make structs of the data
    # ================================================
    # Naming convention: *typeofdata*_*obtainedfrom*_*coordinatesystem* 
    # Example: B_GSM is magnetic field data, from IGRF model, in GSM coordinates
    B_GSM = Struct(B_GSM_long.T,Epoch,None,'nT',['Bx_IGRF (from Earth to Sun)','By IGRF (in magnetic equatorial plane)',"Bz_IGRF (aligned with Earth's mag. dipole axis)",'|B_IGRF|'],'IGRF B-field in GSM')
    B_SC = Struct(B_SC_long.T,Epoch,None,'nT',['Bx (away from Sun)','By','Bz','|B|'],'IGRF B-field in SC coords.')
    E_GSM = Struct(E_GSM_long.T,Epoch,None,'mV/m',['Ex_IGRF','Ey_IGRF','Ez_IGRF','|E|'],'IGRF E-field in GSM')
    E_SC = Struct(E_SC_long.T,Epoch,None,'mV/m',['Ex_IGRF (away from Sun)','Ey_IGRF','Ez_IGRF','|E_IGRF|'],'IGRF E-field in SC coords.')
    v_SC = Struct(v_SC_long.T,Epoch,None,'km/s',['x (away from Sun)','y','z'],'Juice velocity in SC coords.')


    return B_GSM, B_SC, E_GSM, E_SC, v_SC