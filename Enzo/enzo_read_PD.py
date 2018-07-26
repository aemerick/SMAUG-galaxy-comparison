import numpy as np
import matplotlib.pyplot as plt


all_sim_names = ['enzo_equillibrium_noFB_100Mres',
                 'enzo_equillibrium_noFB_1000Mres',
                 'enzo_adiabatic_noFB_100Mres']

#
# Data files are listed as sim name + a Myr designation (0000, 0500, 1000)
#
# Data files contain the phase diagram data (128,128). Field is units of Msun
#
enzo_PD_data = {}
for name in all_sim_names:
    enzo_PD_data[name] = {}

    for time_string in ['0000','0500','1000','1500']:
        enzo_PD_data[name][time_string] = np.loadtxt(name + '_' + time_string + '_PD.dat')


#
# Useage Examples:
#

if __name__ == '__main__':
    extent = [-6,4,0,8]
    vmin   = -2 # log Msun
    vmax   =  7 # log Msun

    # useage example
    #  data = enzo_PD_data['enzo_equillibrium_noFB_1000Mres']['1000']
    #
    #  plt.imshow( np.log10(data.T), extent = extent, interpolation=None, cmap='viridis', origin='lower')
    #  plt.colorbar()
    #  plt.show()


    # plot_everything_example
    xy = (0,7)
    for name in enzo_PD_data.keys():
        times = enzo_PD_data[name].keys()

        fig, ax = plt.subplots(2,2)
        fig.set_size_inches(8,8)
        for axis, t_name in zip( ax.flatten(), np.sort(times)):
            data = enzo_PD_data[name][t_name]
            axis.imshow( np.log10(data.T), extent = extent, interpolation = None, cmap = 'viridis',
                         origin = 'lower', vmin = vmin, vmax = vmax)
            axis.set_xlabel(r'Number Density (cm$^{-3}$)')
            axis.set_ylabel(r'Temperature (K)')
            axis.minorticks_on()
            axis.annotate(t_name + ' Myr', xy=xy, xytext=xy)
        fig.savefig(name + '_PD_panel.png')
        plt.tight_layout()
        plt.close()
#-----
