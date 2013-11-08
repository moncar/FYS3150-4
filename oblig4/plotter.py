# Plotter module
# Reads filenames holding data as argv

from numpy import loadtxt
import matplotlib.pylab as plt
import sys

# Data from files
try:
    filename1 = sys.argv[1]
#    filename2 = sys.argv[2]
#    filename3 = sys.argv[3]
#    filename4 = sys.argv[4]
except IndexError:
    filename1 = 'res.txt'
#    filename2 = 'resfinal2_vec.txt'
#    filename3 = 'resfinal3_vec.txt'
#    filename4 = 'resfinal4_vec.txt'

data1 = loadtxt(filename1);
#data2 = loadtxt(filename2);
#data3 = loadtxt(filename3);
#data4 = loadtxt(filename4);

# Find no. of iterations
N1 = len(data1[:,0])
#N2 = len(data2[:,0])
#N3 = len(data3[:,0])

# Find no. of planets
nK = len(data1[0,:] - 1)/2

print "Length of arrays to plot: ", N1
print "Number of planets to plot: ", nK

if N1 > 1e3:
    slice = 100
else:
    slice = 1


# Find max. error value
# data[x  true u(x)   est v(x)   error]
#for data in [data1, data2, data3]:
#    print("Max. error for N=%6g is:%g" % (len(data[:,0]), max(data[:,3])))
#    print("Min. error for N=%6g is:%g" % (len(data[:,0]), min(data[:,3])))



fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rc({'xtick.labelsize': 20})
#plt.rc({'ytick.labelsize': 20})

ax1 = fig.add_subplot(111)
#ax2 = ax1.twinx()
ax1labels = []
linestyles = ['k-', 'y-', 'b-', 'kx', 'r-', 'y-', 'c-', 'b-', 'g-']

for i in xrange(0,nK-1):
    # Plot planets, last item is Sun
    ax1.plot(data1[::slice,2*i+1], data1[::slice, 2*i+2], linestyles[i])

p19, = ax1.plot(data1[::slice,2*nK-1],data1[::slice,2*nK], 'y-x')  # SUN


ax1.set_xlabel(r'Rel. $x$-position ' + r'$\chi_x$', fontsize=14)
ax1.set_ylabel(r'Rel. $y$-position ' + r'$\chi_y$', fontsize=14)
#ax1.set_ylim([-4, 4])
#ax1.set_xlim([-4, 4])

#ax2.set_ylim([-3, 5])
#ax1grid()
#ax2.grid()
#plt.legend([p11], \
#       [(r'Relative position $\chi$ [AU]') ] )
#        (r'Ground state $n=0 ,\, {\omega}_r =0.5$ '), \
#        (r'Ground state $n=0 ,\, {\omega}_r =1$ '), \
#        (r'Ground state $n=0 ,\, {\omega}_r =5$ ') ] , loc=4 )
#        ('Analytical solution'), \
#        ('Rel. error, N=%g points' % (N1,)),\
#        ('Rel. error, N=%g points' % (N2,)),\
#        ('Rel. error, N=%g points' % (N3,))] )
plt.show()
