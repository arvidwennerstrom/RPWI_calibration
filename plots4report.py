"""
Not neccessary for the rest of the code. This is just a way to easily organize the plots for my report
-Arvid
"""


from E_calibration_main import E_LP_SC, LP_potentials
# Not really part of the script. The plots I used in my thesis report in an ordered way.
# ========================================================================
# plotting_4_report(LP_diffs,LP_diffs_off,LP_potentials,LP_potentials_off,E_LP_SC,E_LP_off_SC,E_IGRF_SC,E_JMAG_SC)





import matplotlib.pyplot as plt

# def plotting_4_report(LP_diffs,LP_diffs_off,LP_potentials,LP_potentials_off,E_LP_SC,E_LP_off_SC,E_IGRF_SC,E_JMAG_SC):
    
plt.figure()
LP_potentials.plot(plt.gca(),None,None,4,10)


plt.show()
