# Plotter module
# Reads filenames holding data as argv

from numpy import loadtxt, zeros, mean
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
N1 = len(data1[0,:])
#N2 = len(data2[:,0])
#N3 = len(data3[:,0])

# Find no. of columns
nK = len(data1[:,0])-1

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
#figE= plt.figure()

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rc({'xtick.labelsize': 20})
#plt.rc({'ytick.labelsize': 20})

ax1 = fig.add_subplot(111)
#axE = figE.add_subplot(111)
#ax2 = ax1.twinx()
ax1labels = []
linestyles = ['k-', 'y-', 'b-', 'kx', 'r-', 'y-', 'c-', 'b-', 'g-']
linestylesX = ['kx', 'yx', 'bx', 'kx', 'rx', 'yx', 'cx', 'bx', 'gx']

for i in xrange(1, nK+1):
#    ax1.plot(data1[0,:], linestylesX[i])
#    ax1.hist(data1[0,:])
    ax1.plot(data1[0,:], data1[i,:]/max(data1[i,:]), linestyles[i+1])
#    ax1labels.append('t = %.3g' % (i*1./(nK-1)))
    #ax1.plot(data1[2*i+1,::slice], data1[2*i+2,::slice], linestyles[i])
    print("Mean: %g" % (mean(data1[i,:]),))

# # Plot errors
#for i in xrange(2, nK-1,  (nK-1)/((nK-1)/42)):
#    err = 1 - abs(  (data1[i,:] + (1 - data2[0,:]) )/(1 - data2[0,:]) )
#    axE.plot(data1[0,:], err) 


ax1.set_xlabel(r'$x$-position ', fontsize=14)
ax1.set_ylabel(r'$u(x,t) = v(x,t) - u_s (x)$ ', fontsize=14)
#axE.set_xlabel(r'$x$-position ', fontsize=14)
#axE.set_ylabel(r'Error $ 1 - |u(x,t) / - u_s (x) |$ ', fontsize=14)
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

box = ax1.get_position()
ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax1.legend(ax1labels, loc='center left', bbox_to_anchor=(1, 0.5))
#axE.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#axE.legend(ax1labels, loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()
