
# Position of Langmuir probes in spacecraft cooridnates [m].
# Rows for probe number [1-4] and columns for SC-axis [X Y Z].
# Source: "The Radio & Plasma Wave Investigation (RPWI)
# for the JUpiter ICy moons Explorer (JUICE), page 33)"

import numpy as np

class RPWI_Data:
    def __init__(self):
        self.LP1_position = 1e-3*np.array([2484, 2895, 4377])
        self.LP2_position = 1e-3*np.array([1455, -3238, 5303])
        self.LP3_position = 1e-3*np.array([1455, -3238, -1849])
        self.LP4_position = 1e-3*np.array([-2768, 2686, 4432])
        self.LP12_distance = self.LP1_position-self.LP2_position
        self.LP13_distance = self.LP1_position-self.LP3_position
        self.LP14_distance = self.LP1_position-self.LP4_position
        self.LP23_distance = self.LP2_position-self.LP3_position
        self.LP24_distance = self.LP2_position-self.LP4_position
        self.LP34_distance = self.LP3_position-self.LP4_position
rpwi_data = RPWI_Data()

calibration_coefficients = np.array([[5.15*1e-6], [4.97*1e-6], [5.07*1e-6], [9.94*1e-5]])