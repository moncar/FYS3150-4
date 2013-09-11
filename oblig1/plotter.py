# Plotter module
# Reads filenames holding data as argv

from numpy import loadtxt
import matplotlib.pylab as plt
import sys

# Data from files
try:
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
except IndexError:
    filename1 = 'res150.txt'
    filename2 = 'res1500.txt'

data1 = loadtxt(filename1);
data2 = loadtxt(filename2);

# Find no. of iterations
N1 = len(data1[:,0])
N2 = len(data2[:,0])

fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rc({'xtick.labelsize': 20})
#plt.rc({'ytick.labelsize': 20})

ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()

p1, = ax1.plot(data1[:,0],data1[:,2], label=('N=%g points' % (N1,)) )
p2, = ax1.plot(data2[:,0],data2[:,2], label=('N=%g points' % (N2,)) )
p3, = ax1.plot(data2[:,0],data2[:,1], label=('Analytical solution') )
p4, = ax2.plot(data1[:,0],data1[:,3], \
   label=('$\log_{10}$ (Rel. error), N=%g points' % (N1,)) ) 
p5, = ax2.plot(data2[:,0],data2[:,3], 
   label=('$\log_{10}$ (Rel. error), N=%g points' % (N2,)) ) 
plt.legend([p1,p2,p3,p4,p5], \
        [('N=%g points' % (N1,)), \
        ('N=%g points' % (N2,)), \
        ('Analytical solution'), \
        ('$\log_{10}$ (Rel. error), N=%g points' % (N1,)),\
        ('$\log_{10}$ (Rel. error), N=%g points' % (N2,))] )
plt.show()
