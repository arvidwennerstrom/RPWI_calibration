import datetime, spiceypy as spice, matplotlib.pyplot as plt


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