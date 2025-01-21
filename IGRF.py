
def IGRF_B_field(rootDir,Epoch):
    import numpy as np
    from geopack import geopack
    from E_calibr_support import Struct, rotation_matrix, tt2000_to_unix
    from load_Ephemeris import load_Ephemeris
    # IGRF B-field using geopack
    # Documentation on geopack: (https://github.com/tsssss/geopack/?tab=readme-ov-file)


    Epoch_ephm, EARTH_SC, pos_ephem_GSM, v_ephem_GSM = load_Ephemeris(rootDir,Epoch)


    # List of rotation matrices.
    R_GSM2SC = np.zeros((3,3,len(Epoch_ephm)))


    # Naming convention: *typeofdata*_*obtainedfrom*_*coordinatesystem* 
    # Example: B_GSM is magnetic field data, from IGRF model, in GSM coordinates
    B_GSM = np.zeros((3,len(Epoch_ephm)))
    v_GSM = v_ephem_GSM


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
        B_GSM[:,i] = bgsm


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
        R_GSM2SC[:,:,i] = R


    # Get E-field in GSM
    # ================================================
    E_GSM = np.cross(-v_GSM.T,B_GSM.T).T



    # Convert to SC coordinates
    # ================================================
    B_SC = np.matmul(R_GSM2SC.transpose(2,0,1),B_GSM[None,...].T)[...,0].T
    E_SC = np.matmul(R_GSM2SC.transpose(2,0,1),E_GSM[None,...].T)[...,0].T
    v_SC = np.matmul(R_GSM2SC.transpose(2,0,1),v_GSM[None,...].T)[...,0].T



    # Add the norm of x,y and z components
    # ================================================
    B_GSM = np.concatenate((B_GSM,np.linalg.norm(B_GSM,axis=0)[None,...]))
    B_SC = np.concatenate((B_SC,np.linalg.norm(B_SC,axis=0)[None,...]))
    E_GSM = np.concatenate((E_GSM,np.linalg.norm(E_GSM,axis=0)[None,...]))
    E_SC = np.concatenate((E_SC,np.linalg.norm(E_SC,axis=0)[None,...]))
    v_GSM = np.concatenate((v_GSM,np.linalg.norm(v_GSM,axis=0)[None,...]))
    v_SC = np.concatenate((v_SC,np.linalg.norm(v_SC,axis=0)[None,...]))


    # Make structs of the data
    # ================================================
    # Naming convention: *typeofdata*_*obtainedfrom*_*coordinatesystem* 
    # Example: B_GSM is magnetic field data, from IGRF model, in GSM coordinates
    B_GSM = Struct(B_GSM,Epoch_ephm,None,'nT',['Bx_IGRF (from Earth to Sun)','By IGRF (in magnetic equatorial plane)',"Bz_IGRF (aligned with Earth's mag. dipole axis)",'|B_IGRF|'],'IGRF B-field in GSM')
    B_SC = Struct(B_SC,Epoch_ephm,None,'nT',['Bx (away from Sun)','By','Bz','|B|'],'IGRF B-field in SC coords.')
    E_GSM = Struct(E_GSM,Epoch_ephm,None,'mV/m',['Ex_IGRF','Ey_IGRF','Ez_IGRF','|E|'],'IGRF E-field in GSM')
    E_SC = Struct(E_SC,Epoch_ephm,None,'mV/m',['Ex_IGRF (away from Sun)','Ey_IGRF','Ez_IGRF','|E_IGRF|'],'IGRF E-field in SC coords.')
    v_SC = Struct(v_SC,Epoch_ephm,None,'km/s',['x (away from Sun)','y','z'],'Juice velocity in SC coords.')


    return B_GSM, B_SC, E_GSM, E_SC, v_SC


    # # Extend length of data to match resolution of LP measurements
    # # ================================================
    # B_GSM_long = np.zeros((len(Epoch),4))
    # B_SC_long = np.zeros((len(Epoch),4))
    # E_GSM_long = np.zeros((len(Epoch),4))
    # E_SC_long = np.zeros((len(Epoch),4))
    # v_SC_long = np.zeros((len(Epoch),4))


    # for i in range(3):
    #     B_GSM_long[:,i] = np.interp(Epoch,Epoch_ephm,B_GSM[:,i])
    #     B_SC_long[:,i] = np.interp(Epoch,Epoch_ephm,B_SC[:,i])
    #     E_GSM_long[:,i] = np.interp(Epoch,Epoch_ephm,E_GSM[:,i])
    #     E_SC_long[:,i] = np.interp(Epoch,Epoch_ephm,E_SC[:,i])
    #     v_SC_long[:,i] = np.interp(Epoch,Epoch_ephm,v_SC[:,i])

