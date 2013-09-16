# Plotter module
# Reads filenames holding data as argv

from numpy import loadtxt
import matplotlib.pylab as plt
import sys

# Data from files
try:
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    filename3 = sys.argv[3]
except IndexError:
    filename1 = 'res10.txt'
    filename2 = 'res100.txt'
    filename3 = 'res1000.txt'

data1 = loadtxt(filename1);
data2 = loadtxt(filename2);
data3 = loadtxt(filename3);

# Find no. of iterations
N1 = len(data1[:,0])
N2 = len(data2[:,0])
N3 = len(data3[:,0])

# Find max. error value
# data[x  true u(x)   est v(x)   error]
for data in [data1, data2, data3]:
    print("Max. error for N=%6g is:%g" % (len(data[:,0]), max(data[:,3])))
    print("Min. error for N=%6g is:%g" % (len(data[:,0]), min(data[:,3])))



fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rc({'xtick.labelsize': 20})
#plt.rc({'ytick.labelsize': 20})

ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()

p11, = ax1.plot(data1[:,0],data1[:,2])
p12, = ax1.plot(data2[:,0],data2[:,2])
p13, = ax1.plot(data3[:,0],data3[:,2])
p2, = ax1.plot(data2[:,0],data2[:,1], color='k', linestyle='-.', linewidth=3)
p31, = ax2.plot(data1[:,0],data1[:,3], linestyle='--')
p32, = ax2.plot(data2[:,0],data2[:,3], linestyle='--') 
p33, = ax2.plot(data3[:,0],data3[:,3], linestyle='--') 
ax1.set_xlabel('x')
ax1.set_ylabel('Estimated derivative')
ax2.set_ylabel(r'\log_{10}$(Rel. error)')
ax1.set_ylim([0, 2])
ax2.set_ylim([-3, 5])
ax1.grid()
ax2.grid()
plt.legend([p11,p12,p13,p2,p31,p32,p33], \
        [('N=%g points' % (N1,)), \
        ('N=%g points' % (N2,)), \
        ('N=%g points' % (N3,)), \
        ('Analytical solution'), \
        ('Rel. error, N=%g points' % (N1,)),\
        ('Rel. error, N=%g points' % (N2,)),\
        ('Rel. error, N=%g points' % (N3,))] )
plt.show()
