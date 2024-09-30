# HK-plotting script
# By: Arvid Wennerström"



# ================================================
#               Imports & Stuff 
# ================================================
import pycdfpp
import os
import re
import numpy as np


# ================================================
#               Main
# ================================================

# Input overall directory path, where folders 
# "datasets", "data_created" and "spice" exist
# ================================================ 
rootDir = "C:/Users/arvwe/Onedrive - KTH/MEX/IRF" # My laptop
rootDir = "C:/Users/arvidwen/Onedrive - KTH/MEX/IRF" # KTH computers



#           Choose data to create
# ================================================ 
date = '240822'
create_LP = True
create_HK10002 = False
create_HK10064 = False



# Create directions and filenames to look for.
# ================================================
datasetDir = rootDir + '/datasets/20' + date[0:2] + '/' + date[2:4] + '/' + date[4:6]
regexPattern = {
    "HK10002":re.compile("^JUICE_LU_RPWI-PPTD-LWYHK10002.*"),
    "HK10064":re.compile("^JUICE_LU_RPWI-PPTD-LWYHK10064.*"),
    "LP":re.compile("^JUICE_L1a_RPWI-LP-SID1_RICH_DE763_SNAP.*")}



# Search directory for files using regular expression and add them to the correct
# "allFiles" dict key.
# ================================================
allFiles = {"HK10002":[],"HK10064":[],"LP":[]}

for filename in os.listdir(datasetDir):
    if regexPattern["HK10002"].match(filename):
        allFiles["HK10002"].append(filename)
    elif regexPattern["HK10064"].match(filename):
        allFiles["HK10064"].append(filename)
    elif regexPattern["LP"].match(filename):
        allFiles["LP"].append(filename)



#       Load and manage LP data
# ================================================
# ================================================
if create_LP:
    mask_description = {
        'explanation': "Mask array is of same size as LP_data, with each index corresponding to a data point. The mask value contains info on data quality, with high values being better quality.",
        '1': "TM saturation (either high or low).",
        '2': "Interference with RIME instrument.",
        '3': "Lunar wake data, will overwrite '1'.",
        '4': "High quality data, this does not contain any contaminations."
        }

    # Max and min TM values --> saturation
    saturation_max_TM = 1048575
    saturation_min_TM = -1048576
    # saturation_max_TM = 1048500
    # saturation_min_TM = -1048500


    file_i = 0
    for filename in allFiles["LP"]:
        #   Load the necessary data from the correct cdf
        # ================================================
        cdf = pycdfpp.load(datasetDir + '/' + filename)
        Epoch_LP = cdf["Epoch"].values['value']
        LP_data = (cdf["LP_DATA"].values).T
        Mask = np.zeros(LP_data.shape)


        # Removes the data from when instrument was in density mode
        # Only applicable for 2024-01-26
        # ================================================
        if True:         
            #Seq. 14: in *T101713 @:            7.59538487e17 < t < 7.5953931e17 
            #Seq. 17: in *T101713 & *T110239 @: 7.5954003e17 < t < 7.59541e17
            #Seq. 18 & 20: in *T110239 @:       t > 7.5954327e17

            # NOTE: The way this is done might be inefficient, with getting a mask and values
            # for all data that is to be removed. It may be faster to check when density mode
            # starts and stops an remove the indeces in between.
            density_mask = np.zeros(len(Epoch_LP))
            if filename == "JUICE_L1a_RPWI-LP-SID1_RICH_DE763_SNAP_20240126T101713_V03.cdf":
                i = 0
                for t in Epoch_LP:
                    if t > 7.59538487*1e17:
                        density_mask[i] = 1      
                        if t < 7.5953931*1e17:
                            density_mask[i] = 1
                        elif t > 7.5954003*1e17:
                            density_mask[i] = 1
                    i+=1
            elif filename == "JUICE_L1a_RPWI-LP-SID1_RICH_DE763_SNAP_20240126T112039_V02.cdf":
                i = 0
                for t in Epoch_LP:
                    if t < 7.59541*1e17:
                        density_mask[i] = 1
                    elif t > 7.5954327*1e17:
                        density_mask[i] = 1
                    i+=1
            Epoch_LP = Epoch_LP[density_mask == 0]
            LP_data = LP_data[:,density_mask == 0]
            Mask = Mask[:,density_mask == 0]
       


        # Removes config. start noise from certain data files
        # ================================================
        if True:
            list_of_start_files = [
                "JUICE_L1a_RPWI-LP-SID1_RICH_DE763_SNAP_20240126T090524_V03",
                "JUICE_L1a_RPWI-LP-SID1_RICH_DE763_SNAP_20240706T023204_V01.cdf",
                "JUICE_L1a_RPWI-LP-SID1_RICH_DE763_SNAP_20240819T202634_V01.cdf",
                "JUICE_L1a_RPWI-LP-SID1_RICH_DE763_SNAP_20240820T180737_V01.cdf",
                "JUICE_L1a_RPWI-LP-SID1_RICH_DE763_SNAP_20240823T035821_V01.cdf",
                "JUICE_L1a_RPWI-LP-SID1_RICH_DE763_SNAP_20240126T090524_V03.cdf"]
            
            for filename_2 in list_of_start_files:    
                if filename == filename_2:
                    removal_width = 32 # Number of samples to remove
                    Epoch_LP = Epoch_LP[removal_width-1:]
                    LP_data = LP_data[:,removal_width-1:]
                    Mask = Mask[:,removal_width-1:]


    
        # Fills masking with saturations
        # ================================================
        if True:
            # When LP_data reaches TM saturation values (max or min) for any channel,
            # set corresponding mask value to 2
            Mask[LP_data >= saturation_max_TM] = 1
            Mask[LP_data <= saturation_min_TM] = 1



        # Add interference (RIME etc.) to mask
        # ================================================
        if True:
            # NOTE: Work in progress
            pass



        # Add Lunar wake data to mask 
        # Only applicable on 2024-08-19
        # ================================================
        if True: 
            # This file contains lunar wake data
            if filename == "JUICE_L1a_RPWI-LP-SID1_RICH_DE763_SNAP_20240819T202634_V01.cdf":
                i = 0
                for t in Epoch_LP:
                    if t > 7.77371896287*1e17:   
                        if t < 7.77373767102*1e17:
                            Mask[:,i] = 3
                    i+=1



        # Data that has not been filtered out gets the 
        # highest flag value (most valuable data)
        # ================================================
        Mask[Mask == 0] = 4

    
        # Combine previous and current files.
        # ================================================
        if file_i == 0:
            Epoch_LP_All = Epoch_LP
            LP_data_All = LP_data
            Mask_All = Mask
        else:
            Epoch_LP_All = np.concatenate((Epoch_LP_All,Epoch_LP))
            LP_data_All = np.concatenate((LP_data_All,LP_data),axis = 1)
            Mask_All = np.concatenate((Mask_All,Mask),axis = 1)


        file_i +=1
        #           End of loop over LP files
        # ================================================


    

    # Removes config. change noise from all data
    # ================================================
    if True:
        removal_width = 32 # Number of samples to remove
        new_config_delay = 1e8 # [ns]  new_config_delay at which a new config. is assumed to have been made

        noise_mask = np.zeros(len(Epoch_LP_All))
        # Mask the next samples when dt is above new_config_delay, for which a new config. has been made
        delta_t = Epoch_LP_All[1:] - Epoch_LP_All[:-1]

        i = 0
        for dt in delta_t:
            if dt > new_config_delay:
                noise_mask[i:i+removal_width] = 1
            i +=1

        Epoch_LP_All = Epoch_LP_All[noise_mask == 0]
        LP_data_All = LP_data_All[:,noise_mask == 0]
        Mask_All = Mask_All[:,noise_mask == 0]


    #           Save the data in a .npz file
    # ================================================
    np.savez((rootDir + "/data_created/LP-SID1_" + date), Epoch=Epoch_LP_All, LP_data=LP_data_All, Mask = Mask_All)



#       Handling HK10002 data, LP preamp. temp. 
# ================================================
# ================================================
if create_HK10002:
    file_i = 0
    for file in allFiles["HK10002"]:

        #   Load the necessary data from the correct cdf
        # ================================================
        cdf = pycdfpp.load(datasetDir + '/' + file)
        Epoch_HK10002 = cdf["Epoch"].values['value']
        LWT0449D = cdf["LWT0449D_CALIBRATED"].values
        LWT0449E = cdf["LWT0449D_CALIBRATED"].values
        LWT0449F = cdf["LWT0449D_CALIBRATED"].values
        LWT044A0 = cdf["LWT0449D_CALIBRATED"].values
        HK10002_data = [LWT0449D, LWT0449E, LWT0449F, LWT044A0]
        

        # Combine previous and current files.
        # ================================================
        if file_i == 0:
            Epoch_HK10002_All = Epoch_HK10002
            HK10002_data_All = HK10002_data
        else:
            Epoch_HK10002_All = np.concatenate((Epoch_HK10002_All,Epoch_HK10002))
            HK10002_data_All = np.concatenate((HK10002_data_All,HK10002_data),axis = 1)
            

        file_i +=1
        #           End of loop over HK10002 files
        # ================================================

    #           Save the data in a .npz file
    # ================================================
    np.savez((rootDir + "/data_created/LWYHK10002_" + date),Epoch=Epoch_HK10002_All, HK10002_data=Epoch_HK10002_All)



#       Handling HK10064 data, E-box temp. 
# ================================================
# ================================================
if create_HK10064:
    file_i = 0
    for file in allFiles["HK10064"]:

        #   Load the necessary data from the correct cdf
        # ================================================
        cdf = pycdfpp.load(datasetDir + '/' + file)
        Epoch_HK10064 = cdf["Epoch"].values['value']
        HK10064_data = cdf["LWT04567"].values
        

        # Combine previous and current files.
        # ================================================
        if file_i == 0:
            Epoch_HK10064_All = Epoch_HK10064
            HK10064_data_All = HK10064_data
        else:
            Epoch_HK10064_All = np.concatenate((Epoch_HK10064_All,Epoch_HK10064))
            HK10064_data_All = np.concatenate((HK10064_data_All,HK10064_data))
            

        file_i +=1
        #           End of loop over HK10064 files
        # ================================================

    #           Save the data in a .npz file
    # ================================================
    np.savez((rootDir + "/data_created/LWYHK10064_" + date),Epoch=Epoch_HK10064_All, HK10064_data=Epoch_HK10064_All)
