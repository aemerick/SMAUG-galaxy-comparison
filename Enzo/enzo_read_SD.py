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
enzo_SD_data = {}
for name in all_sim_names:
    enzo_SD_data[name] = {}

    for axis_name in ['x','z']:
        enzo_SD_data[name][axis_name] = {}

        for time_string in ['0000','0500','1000','1500']:
            enzo_SD_data[name][axis_name][time_string] =\
                      np.loadtxt(name + '_' + time_string + '_' + axis_name + '_SD.dat')


#
# Useage Examples:
#

if __name__ == '__main__':
    extent = [-5,5,-5,5]
    vmin   = -1.5 # log Msun
    vmax   =  1.5 # log Msun

    # useage example
    #  data = enzo_SD_data['enzo_equillibrium_noFB_1000Mres']['1000']
    #
    #  plt.imshow( np.log10(data), extent = extent, interpolation=None, cmap='viridis', origin='lower')
    #  plt.colorbar()
    #  plt.show()


    # plot_everything_example
    for name in enzo_SD_data.keys():
        for axname in enzo_SD_data[name].keys():
            times = enzo_SD_data[name][axname].keys()

            fig, ax = plt.subplots(2,2, sharex=True)
            fig.set_size_inches(8,8)
            for axis, t_name in zip( ax.flatten(), np.sort(times)):
                data = enzo_SD_data[name][axname][t_name]
                axis.imshow( np.log10(data), extent = extent, interpolation = None, cmap = 'viridis',
                             origin = 'lower', vmin = vmin, vmax = vmax)
#            axis.set_xlabel(r'Number Density (cm$^{-3}$)')
#            axis.set_ylabel(r'Temperature (K)')
#            axis.minorticks_on()
#            axis.annotate(t_name + ' Myr', xy=xy, xytext=xy)
                axis.set_xticklabels([])
                axis.set_yticklabels([])
            fig.subplots_adjust(hspace = 0, wspace = 0)
            fig.savefig(name + '_' + axname + '_SD_panel.png')
            plt.tight_layout()
            plt.close()
#-----
