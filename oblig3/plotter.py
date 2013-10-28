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

print "Length of arrays to plot: ", N1

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

p10, = ax1.plot(data1[::slice,1], data1[::slice, 2], 'k-')
p11, = ax1.plot(data1[::slice,3], data1[::slice, 4], 'y-')
p12, = ax1.plot(data1[::slice,5], data1[::slice, 6], 'b-')
#p13, = ax1.plot(data1[::slice,7], data1[::slice, 8], 'kx')
p14, = ax1.plot(data1[::slice,9], data1[::slice,10], 'r-')
#p13, = ax1.plot(data1[::slice,7], data1[::slice, 8], 'r-')
#p14, = ax1.plot(data1[::slice,9], data1[::slice,10], 'c-')
p15, = ax1.plot(data1[::slice,11],data1[::slice,12], 'y-')
p16, = ax1.plot(data1[::slice,13],data1[::slice,14], 'g-')
p17, = ax1.plot(data1[::slice,15],data1[::slice,16], 'b-')
p18, = ax1.plot(data1[::slice,17],data1[::slice,18], 'g-')
#p19, = ax1.plot(data1[::slice,17],data1[::slice,18], 'y-x')  # SUN
p19, = ax1.plot(data1[::slice,19],data1[::slice,20], 'y-x')  # SUN

#p2, = ax1.plot(data2[:,0],data2[:,1], color='k', linestyle='-.', linewidth=3)
#p31, = ax2.plot(data1[:,0],data1[:,1], linestyle='--')
#p32, = ax2.plot(data2[:,0],data2[:,1], linestyle='--') 
#p33, = ax2.plot(data3[:,0],data3[:,1], linestyle='--') 
ax1.set_xlabel(r'Rel. $x$-position ' + r'$\chi_x$', fontsize=14)
ax1.set_ylabel(r'Rel. $y$-position ' + r'$\chi_y$', fontsize=14)
#ax2.set_ylabel(r'\log_{10}$(Rel. error)', fontsize=14)
#ax1.set_ylim([0, 2])
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
