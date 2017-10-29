#!/usr/bin/env python
"""
This program plots the Hagedron density tau for given B,S,I values.
"""

from __future__ import print_function
import argparse

import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
#from scipy import optimize


import ReadArr

print(ReadArr.__doc__)
print(ReadArr.read.__doc__)

ReadArr.read.readarrays()

def doPlot( iB, iS, iI ):
    """Do the plot
    """

    print("BSI:",iB,iS,iI)


    if ((abs(iB)<2)&(abs(iS)<4)&(iI<4)):
        xvals1 = np.arange(1,300+1)*ReadArr.read.dm
        yvals1 = ReadArr.read.arrtauhadron[:,iB,iS,iI]
        #    print(yvals1)
        line, = plt.plot(xvals1,yvals1,"r:")
#        line.set_antialiased(False)

    xvals2 = np.arange(1,ReadArr.read.nbin+1)*ReadArr.read.dm
    yvals2 = ReadArr.read.arrtau[:,iB,iS,iI]

    line, = plt.plot(xvals2,yvals2,"b--")
#    line.set_antialiased(False)
#    print(yvals2)

    line, = plt.plot(xvals2,np.exp(xvals2/0.16),"g:")

    ax = plt.gca()
    plt.yscale("log")
    plt.ylim(ymin=1e-2)
    plt.xlabel('M [GeV]')
    plt.text(0.7,0.1,'B=%d, S=%#d, I=%2.1f' % (iB,iS,iI/2.),
             horizontalalignment='center',
             verticalalignment='center',
             transform=ax.transAxes,
             fontsize=28)

    plt.savefig('fig%+02d%+02d+y%02d.eps'%(iB,iS,iI))

    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot Tau.',
        epilog=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--qB','-b',
                        type=int, default=0,
                        help='baryon number (default: %(default)s)')
#    parser.add_argument('--qQ','-q',
#                        type=int, default=0,
#                        help='charge (default: %(default)s)')
    parser.add_argument('--qS','-s',
                        type=int, default=0,
                        help='strangeness (default: %(default)s)')
    parser.add_argument('--qI','-i',
                        type=int, default=2,
                        help='isospin (default: %(default)s)')

    args = parser.parse_args()

    doPlot(args.qB,
           args.qS,
           args.qI)
