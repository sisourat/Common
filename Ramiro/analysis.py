import sys
import numpy as np
import matplotlib.pyplot as plt


dat = np.loadtxt(sys.argv[1])
fsta = dat[:,0]
esta = dat[:,1]

prob = float(np.count_nonzero(fsta == 2))/len(fsta)

#hist.append(wigner(2.0,0.0,0.28)[1])
plt.hist(esta, bins='auto')  # arguments are passed to np.histogram
plt.show()

print(fsta,prob)
