import numpy as np
import glob
import sys

def extract_data(name, output='vcirc.dat'):
    """
    Read MakeDiskGalaxy output - move to two column
    """

    print "Reading " + name

    data = np.genfromtxt(name)

    npoints = int(data[0])

    r  = data[1:npoints+1]
    v2 = data[npoints + 1: npoints + npoints + 1]
    v  = np.sqrt(v2)

    np.savetxt(output, np.array([r,v]).T, fmt="%8.8E")

    print "Wrote " + output

    return


if __name__ == '__main__':

    if len(sys.argv) == 1:
        velocity_file = glob.glob('*.vc')[0]

        extract_data(velocity_file)

    elif len(sys.argv) == 2:

        extract_data(sys.argv[1])

    elif len(sys.argv) == 3:

        extract_data(sys.argv[1], output = sys.argv[2])
