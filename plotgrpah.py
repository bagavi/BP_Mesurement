import numpy
from scipy import *
from scipy import signal
from matplotlib.pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.mlab as ax
from scipy.linalg import lstsq
import scipy.special as sp

P1 = loadtxt("testdatafrompython.dat")
P2 = loadtxt("testdatafrompythonheart.dat")
B = loadtxt("outputtest.dat")
C = loadtxt("outputtest2.dat")
D = loadtxt("outputtest3.dat")
E = loadtxt("outputtest4.dat")
plot(B ,C )

plot(B[0::2] ,C[0::2] )
plot(B[1::2] ,C[1::2] )
plot(B[0::2] ,E[0::2] )
plot(B[1::2] ,E[1::2] )
figure()
plot(D-8.13)
title("from the c code")
figure()

#title("difference in the values")
len1 = min (len(P1) , len(B))
plot(P2)
title("from the python code")
show()

