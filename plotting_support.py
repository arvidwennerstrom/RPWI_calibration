import datetime, spiceypy as spice, matplotlib.pyplot as plt, numpy as np
from matplotlib.ticker import FuncFormatter, MaxNLocator


# Does not do anything rn, might be useful in the future
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


    def plot(self, axis, data_range = None, chunk = None, mask_limit = 4, plot_style = {'linewidth':'1'}):
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
            [axis.plot(time,data[i], **plot_style, label=self.legends[i]) for i in data_range]
        else:
            [axis.plot(time[mask[i] >= mask_limit],data[i][mask[i] >= mask_limit], **plot_style, label=self.legends[i]) for i in data_range]

        axis.set_xlabel('Timestamp (UTC)')
        axis.set_ylabel(self.unit)
        axis.set_title(self.title)
        axis.legend(loc='upper right')
        axis.grid(True)
              
        # Ensures there are more labeled points on the time axis 
        # and converts them to human-readable timestamps. 
        # Also rotates them for readablility
        axis.xaxis.set_major_locator(MaxNLocator(nbins=10))
        axis.xaxis.set_major_formatter(FuncFormatter(dynamic_time_formatter))
        plt.gcf().autofmt_xdate()

        return axis


def datenum_to_tt2000(datenum):
    # MATLAB's datenum is the number of days since 0-Jan-0000
    # So, adjust by the difference between Python's datetime epoch (1970-01-01) and MATLAB's epoch
    matlab_epoch = datetime.datetime(1970, 1, 1)
    days = datenum - 719529  # MATLAB's epoch is 719529 days before 1970-01-01
    dt = matlab_epoch + datetime.timedelta(days=days)

    utc_str = dt.strftime('%Y-%m-%dT%H:%M:%S.%f')
    
    et = spice.str2et(utc_str)
    # Convert Ephemeris Time to TT2000 (you might need to load necessary kernels)
    tt2000_time = spice.unitim(et, 'ET', 'TT')
    return tt2000_time



def tt2000_to_readable(tt2000,precision = 9):    
    # Convert TT2000 (nanoseconds) to ephemeris time (TDB) in seconds
    ephem_time = spice.unitim(tt2000/ 1e9, 'TT', 'TDB')

    # Convert ephemeris time to human-readable UTC string
    utcs = spice.spiceypy.et2utc(ephem_time , format_str='ISOC', prec=precision)

    readable = utcs
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