
import math, numpy as np

def fit_to_E_obs(Epoch, E_obs, E_vxB, Vs, cross_E_coeffs = False):
    
    """
    Find a function that approximates measured E-field, E_obs. Possible contributing factors include:
        - E_vxB, the expected E-field derived from E = -v x B, where v is plasma velocity w.r.t the SC, and B is magnetic field.
        - Vs, spacecraft potential

    """

    # Step configuration. This
    duration = (Epoch[-1] - Epoch[0])*1e-9
    max_step_duration = 20*60 # [sec]
    n_steps = math.ceil(duration/max_step_duration)
    step_size = math.ceil(len(Epoch)/n_steps)

    overlap_fraction = 0.9
    window_shift = int((1 - overlap_fraction) * step_size)


  # Generate overlapping windows
    overlapping_windows = []
    start_idx = 0
    while start_idx < len(Epoch):
        end_idx = start_idx + step_size-1
        if end_idx > len(Epoch):
            end_idx = len(Epoch)
        overlapping_windows.append((start_idx, end_idx))
        start_idx += window_shift


    # Create list of E-field windows
    E_vxBx = [E_vxB[0][start:end] for start, end in overlapping_windows]
    E_vxBy = [E_vxB[1][start:end] for start, end in overlapping_windows]
    E_vxBz = [E_vxB[2][start:end] for start, end in overlapping_windows]
    Vs = [Vs[start:end] for start, end in overlapping_windows]
    E_obs_slices = [E_obs[:, start:end] for start, end in overlapping_windows]
    list_of_fit_Epochs = [(Epoch[start] + Epoch[end-1]) / 2 for start, end in overlapping_windows]
    

    # Create coefficients. They are of length 4 to correspond to each LP channel: P12, P23, P34 and P4
    a = np.zeros((len(list_of_fit_Epochs),3,3))
    b = np.zeros((len(list_of_fit_Epochs),3))
    c = np.zeros((len(list_of_fit_Epochs),3))


    E_fit_combined = np.zeros_like(E_obs)
    weights = np.zeros_like(E_obs)


    for step, (start, end) in enumerate(overlapping_windows):
        for axis in range(3):

            # Allows cross connection between E-components, eg. x-axis can be affected by values in y and z
            if cross_E_coeffs:
                X = np.vstack([E_vxBx[step], E_vxBy[step], E_vxBz[step], Vs[step], np.ones(len(Vs[step]))]).T
                coeffs, residual, rank, s = np.linalg.lstsq(X, E_obs_slices[step][axis], rcond=None)
                a[step][axis][0], a[step][axis][1], a[step][axis][2], b[step][axis], c[step][axis] = coeffs


                # Compute the fitted E for this step
                E_fit_step = (a[step][axis][0] * E_vxBx[step] +
                      a[step][axis][1] * E_vxBy[step] +
                      a[step][axis][2] * E_vxBz[step] +
                      b[step][axis] * Vs[step] +
                      c[step][axis])
                
            # No cross connection allowed
            else:
                E_vxB = [E_vxBx, E_vxBy, E_vxBz][axis]
                X = np.vstack([E_vxB[step], Vs[step], np.ones(len(Vs[step]))]).T
                coeffs, residual, rank, s = np.linalg.lstsq(X, E_obs_slices[step][axis], rcond=None)
                a[step][axis][axis], b[step][axis], c[step][axis] = coeffs
                

                # Compute the fitted E for this step
                E_fit_step = (a[step][axis][axis] * E_vxB[step] +
                      b[step][axis] * Vs[step] +
                      c[step][axis])
                

            # Add E of this step to the combined E. 
            # NOTE: This is now done uniformly, possibly multiply with Gaussian function
            # or similar to smoothen the transition between steps 
            E_fit_combined[axis, start:end] += E_fit_step
            weights[axis, start:end] += 1


    E_fit_combined = E_fit_combined/weights


    return E_fit_combined, a, b, c, list_of_fit_Epochs



# NOTE: Old version, as a backup. Mainly this has no overlap possibility
"""
def fit_to_E_obs(Epoch, E_obs, E_vxB, Vs, cross_E_coeffs = False):
    

    Find a function that approximates measured E-field, E_obs. Possible contributing factors include:
        - E_vxB, the expected E-field derived from E = -v x B, where v is plasma velocity w.r.t the SC, and B is magnetic field.
        - Vs, spacecraft potential



    # Step configuration
    duration = (Epoch[-1] - Epoch[0])*1e-9
    max_step_duration = 20*60
    n_steps = max(math.floor(duration/max_step_duration),1)


    # Create coefficients. They are of length 4 to correspond to each LP channel: P12, P23, P34 and P4
    a = np.zeros((n_steps,3,3))
    b = np.zeros((n_steps,3))
    c = np.zeros((n_steps,3))

    # These are the median time stamps of the steps. Coeffiecients will correspond to one of these.
    list_of_fit_Epochs = []
    [list_of_fit_Epochs.append((t[-1]+t[0])/2) for t in np.array_split(Epoch,n_steps)]
    list_of_fit_Epochs = np.array(list_of_fit_Epochs)


    E_vxBx = np.array_split(E_vxB[0],n_steps)
    E_vxBy = np.array_split(E_vxB[1],n_steps)
    E_vxBz = np.array_split(E_vxB[2],n_steps)
    Vs = np.array_split(Vs,n_steps)

    for axis in range(3):
        for step in range(n_steps):
            # General process:
            # Make the matrix X, to represent the values of each parameter and solve 
            # X using least square to find optimal coefficeints. Use coefficients 
            # to get an E fit to E_obs, for the current step. These are then combined 
            # togheter to cover the whole axis, and then repeated for all axes x,y,z


            # Allows cross connection between E-components, eg. x-axis can be affected by values in y and z
            if cross_E_coeffs:
                X = np.vstack([E_vxBx[step], E_vxBy[step], E_vxBz[step], Vs[step], np.ones(len(Vs[step]))]).T
                coeffs, residual, rank, s = np.linalg.lstsq(
                    X,np.array_split(E_obs[axis], n_steps)[step], rcond=None)
                
                a[step][axis][0], a[step][axis][1], a[step][axis][2], b[step][axis], c[step][axis] = coeffs            
                E_fit_step = a[step][axis][0]*E_vxBx[step] + a[step][axis][1]*E_vxBy[step] + a[step][axis][2]*E_vxBz[step] + b[step][axis]*Vs[step] + c[step][axis]
                
            # No cross connection allowed
            else:
                X = np.vstack([[E_vxBx, E_vxBy, E_vxBz][axis][step], Vs[step], np.ones(len(Vs[step]))]).T
                coeffs, residual, rank, s = np.linalg.lstsq(
                    X,np.array_split(E_obs[axis], n_steps)[step], rcond=None)
                
                a[step][axis][axis], b[step][axis], c[step][axis] = coeffs
                E_fit_step = a[step][axis][axis]*[E_vxBx, E_vxBy, E_vxBz][axis][step] + b[step][axis]*Vs[step] + c[step][axis]
        
    
            # Combine the steps togheter
            if step == 0:
                E_fit_axis = E_fit_step
            else:
                E_fit_axis = np.concatenate((E_fit_axis,E_fit_step))

        # Stack x,y and z into one np.array
        if axis == 0:
            E_fit = E_fit_axis
        else:
            E_fit = np.column_stack((E_fit,E_fit_axis))
    
    
    return E_fit.T, a, b, c, list_of_fit_Epochs
"""