import numpy as np
import matplotlib.pyplot as plt


#
# Files contain two named columns:
#   1) Left edge of time bins (in Myr) - 10 Myr bins
#   2) SFR in Msun/yr
#
all_filenames = ['enzo_equillibrium_noFB_100Mres_SFR.dat',  # 100  Msun gas refinement
                 'enzo_equillibrium_noFB_1000Mres_SFR.dat'] # 1000 Msun gas refinement


enzo_SFR_data = {}
for filename in all_filenames:
    dict_name = filename.strip('_SFR.dat')
    enzo_SFR_data[dict_name] = np.genfromtxt(filename, names = True)


if __name__ == "__main__":
    #
    # only make this plot if script called from command line
    #
    _temp = enzo_SFR_data['enzo_equillibrium_noFB_100Mres']
    plt.plot(_temp['time'], _temp['SFR'])

    _temp = enzo_SFR_data['enzo_equillibrium_noFB_1000Mres']
    plt.plot(_temp['time'], _temp['SFR'])
    plt.yscale('log')
    plt.show()
