import copy
import math
import cmath
import random
import sys
import datetime
import time
import matplotlib.pyplot as plt
import os

# =========================================================0

np = True
try:
    import numpy
except ImportError:
    print('numpy package not installed')
    np = False

# some constants
DEBUG = False
CM_TO_HARTREE = 1. / 219474.6  # 4.556335252e-6 # conversion factor from cm-1 to Hartree
HARTREE_TO_EV = 27.211396132    # conversion factor from Hartree to eV
U_TO_AMU = 1. / 5.4857990943e-4            # conversion from g/mol to amu
ANG_TO_BOHR = 1. / 0.529177211  # 1.889725989      # conversion from Angstrom to bohr
PI = math.pi

def wigner(Q0, P0, fwhm):
    """This function calculates the Wigner distribution for
       a single one-dimensional harmonic oscillator."""
    fwhm=fwhm/HARTREE_TO_EV
    alp = fwhm/2.355
    while True:
     Q = random.random() * 2.0 + 1.07
     P = random.random() * 1000.0 - 500.0
     wig = math.exp(-(Q-Q0)**2/alp**2) * math.exp(-(P-P0)**2*alp**2)/PI
     if wig>1 or wig<0 :
         pass
     elif wig> random.random():
         break

    return Q, P

def main():
    '''Main routine'''

    hist = []
    for i in range(10):
        f = open("init0", "w")
        print("6077.8",wigner(2.0,0.0,0.28)[0],wigner(2.0,0.0,0.28)[1], file = f)
        #print("6077.8  2.08 0.00", file = f)
        print("1.0 -100.0 1.10", file = f)
        f.close()
        os.system("./main.exe input")
        os.system("cat fort.101 >> prob")
        os.system("cat fort.200 >> traj"+str(i))


#        hist.append(wigner(2.0,0.0,0.28)[1])
#    plt.hist(hist, bins='auto')  # arguments are passed to np.histogram
#    plt.show()


if __name__ == "__main__":
    main()
