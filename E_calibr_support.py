import datetime, spiceypy as spice, matplotlib.pyplot as plt, numpy as np


"""     Data structure
========================================================================
Allows saving data in a way where it can easily be retrieved with correct 
untis and explanation. Also includes a plotting function.

"""
class Struct:
    def __init__(self,data,time,mask,unit,legends,title,explaination = None):
        self.data = data
        self.time = time
        self.mask = mask
        self.unit = unit
        self.legends = legends
        self.title = title
        self.explaination = explaination


    def __getitem__(self,item):
        return self.data[item]


    def plot(self, fig_axis, data_range = None, chunk = None, mask_limit = 4, xtick_in_minutes = 0, plot_style = {'linewidth':'1'}):
        from matplotlib.ticker import FuncFormatter, MaxNLocator
        
        time = self.time
        data = self.data
        mask = self.mask

        # Allows plotting of some of the data
        if data_range is None:
            data_range = range(len(self.data))

        # Select only data as part of this chunk
        if chunk is not None:
            time = time[chunk]
            data = data[:,chunk]
            if mask is not None:        
                mask = mask[:,chunk]    


        # Apply masking if it exists
        if self.mask is None:
            [fig_axis.plot(time,data[i], **plot_style, label=self.legends[i]) for i in data_range]
        else:
            [fig_axis.plot(time[mask[i] >= mask_limit],data[i][mask[i] >= mask_limit], **plot_style, label=self.legends[i]) for i in data_range]

        fig_axis.set_xlabel('Timestamp (UTC)')
        fig_axis.set_ylabel(self.unit)
        fig_axis.set_title(self.title)
        fig_axis.legend(loc='upper right')
        fig_axis.grid(True)


        # If xtick size has been specified, this will display them correctly
        if xtick_in_minutes != 0: 
            start = tt2000_to_readable(time[0])
            step_size_tt2000 = xtick_in_minutes*60*1e9

            xticks = np.arange(time[0], time[-1], step_size_tt2000)
            fig_axis.set_xticks(xticks)
            fig_axis.xaxis.set_major_formatter(FuncFormatter(dynamic_time_formatter))

        else:
            # Ensures there are more labeled points on the time axis 
            # and converts them to human-readable timestamps. 
            fig_axis.xaxis.set_major_locator(MaxNLocator(nbins=10))
            fig_axis.xaxis.set_major_formatter(FuncFormatter(dynamic_time_formatter))
        
        # Rotate them for readablility
        plt.gcf().autofmt_xdate()


        return fig_axis



"""     RPWI Data
========================================================================
"""
# Position of Langmuir probes in spacecraft cooridnates [m].
# Rows for probe number [1-4] and columns for SC-axis [X Y Z].
# Source: "The Radio & Plasma Wave Investigation (RPWI)
# for the JUpiter ICy moons Explorer (JUICE), page 33)"
class RPWI_Data:
    def __init__(self):
        self.LP1_position = 1e-3*np.array([2484, 2895, 4377])   # [m]
        self.LP2_position = 1e-3*np.array([1455, -3238, 5303])  # [m]
        self.LP3_position = 1e-3*np.array([1455, -3238, -1849]) # [m]
        self.LP4_position = 1e-3*np.array([-2768, 2686, 4432])  # [m]
        self.LP12_distance = self.LP1_position-self.LP2_position    # [m]
        self.LP13_distance = self.LP1_position-self.LP3_position    # [m]
        self.LP14_distance = self.LP1_position-self.LP4_position    # [m]
        self.LP23_distance = self.LP2_position-self.LP3_position    # [m]
        self.LP24_distance = self.LP2_position-self.LP4_position    # [m]
        self.LP34_distance = self.LP3_position-self.LP4_position    # [m]

        # M is transformation matrix in: U = M*E, where U is voltage and E is electric field
        self.M_E2U = np.array([self.LP12_distance, self.LP23_distance, self.LP34_distance]) # [V/m] --> [V]
        self.M_U2E = np.linalg.inv(self.M_E2U)  # [V] --> [V/m]

rpwi_data = RPWI_Data()


# Lars' temperature coefficients, to convert from TM units to voltage 
coeffs_TM2voltage = np.array([[5.15*1e-6], [4.97*1e-6], [5.07*1e-6], [9.94*1e-5]])



"""     Time conversions
========================================================================
These are used to convert between different time formats.
"""
def year_doy_to_tt2000(year, doy, hour):
    # Convert year and DOY to a datetime object
    base_date = datetime.datetime(year, 1, 1) + datetime.timedelta(days=doy - 1, hours=hour)
    
    # TT2000 reference point: 2000-01-01T12:00:00
    tt2000_reference = datetime.datetime(2000, 1, 1, 12, 0, 0)
    
    # Calculate the difference in nanoseconds
    delta = base_date - tt2000_reference
    tt2000 = delta.total_seconds() * 1e9  # Convert to nanoseconds
    return int(tt2000)



def tt2000_to_unix(tt2000):
    """
    Convert TT2000 time to Unix time.
    Parameters:
        tt2000 (int): TT2000 time in nanoseconds since 2000-01-01T12:00:00.000000000.
    Returns:
        float: Unix time in seconds since 1970-01-01T00:00:00.
    """
    # Convert TT2000 (nanoseconds) to ephemeris time (ET/TDB) in seconds   
    ephem_time = spice.unitim(tt2000*1e-9, 'TT', 'ET')

    # Seconds between TT2000 and Unix epochs (should be 946727935.0)
    ephem_unix_diff = -spice.str2et('1970-01-01T00:00:00') 
    
    # Adjust to Unix epoch
    unix_time = ephem_time + ephem_unix_diff
    
    return unix_time



def datenum_to_tt2000(datenum):
    # NOTE: NOT IN USE ANYMORE
    #  
    # MATLAB's datenum is the number of days since 0-Jan-0000
    # So, adjust by the difference between Python's datetime epoch (1970-01-01) and MATLAB's epoch
    
    matlab_epoch = datetime.datetime(1970, 1, 1)
    days = datenum - 719529  # MATLAB's epoch is 719529 days before 1970-01-01
    dt = matlab_epoch + datetime.timedelta(days=days)

    utc_str = dt.strftime('%Y-%m-%dT%H:%M:%S.%f')
    
    et = spice.str2et(utc_str)

    # Convert Ephemeris Time to TT2000 (might need to load necessary kernels)
    tt2000_time = spice.unitim(et, 'ET', 'TT')
    return tt2000_time



def tt2000_to_readable(tt2000,precision = 9):    
    # Convert TT2000 (nanoseconds) to ephemeris time (ET/TDB) in seconds   
    ephem_time = spice.unitim(tt2000*1e-9, 'TT', 'ET')

    # Convert ephemeris time to human-readable UTC string
    readable = spice.et2utc(ephem_time , format_str='ISOC', prec=precision)
    return readable



# Custom formatter function to convert displayed tt2000 values dynamically
def dynamic_time_formatter(x, pos):
    # Get the current limits of the x-axis to adjust precision based on zoom level
    current_xlim = plt.gca().get_xlim()
    range_x = current_xlim[1] - current_xlim[0]

    # Adjust precision based on zoom level
    if range_x > 1e12:  # Very zoomed out, show only date
        precision = 0
    elif range_x > 1e9:  # Show hours
        precision = 3
    elif range_x > 1e6:  # Show hours and minutes
        precision = 6
    else:  # Show full precision down to seconds
        precision = 9

    # Convert the tt2000 time to human-readable format with the determined precision
    readable_time = tt2000_to_readable(x,precision)
    return readable_time


"""     Rotation matrix 
========================================================================
"""
def rotation_matrix(v1,v2):
    # Rotation matrix from v1 to v2, by Rodrigues rotation formula

    # Ensure both vectors are unit vectors
    v1 = v1/np.linalg.norm(v1)
    v2 = v2/np.linalg.norm(v2)

    # Compute the axis of rotation (cross product)
    axis = np.cross(v1, v2)
    
    # Compute the angle of rotation (dot product and arccosine)
    cos_theta = np.dot(v1, v2)
    theta = np.arccos(cos_theta)  # Angle in radians
    
    # If the vectors are already aligned (angle is 0), return the identity matrix
    if np.isclose(theta, 0):
        return np.eye(3)
    
    # If the vectors are opposite, we need a special case
    if np.isclose(theta, np.pi):
        # Find an orthogonal vector to use as rotation axis
        axis = np.array([1, 0, 0]) if not np.allclose(v1, [1, 0, 0]) else np.array([0, 1, 0])
        # Use the perpendicular axis to rotate by 180 degrees
        return -np.eye(3) + 2 * np.outer(axis, axis)

    # Normalize the rotation axis
    axis = axis / np.linalg.norm(axis)
    
    # Compute components of Rodrigues' rotation formula
    K = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])

    # Construct the rotation matrix
    R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * np.matmul(K, K)

    return R

