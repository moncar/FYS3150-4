# Moon analyser module
# Reads filename holding data as argv[1]
#       first column holding x pos of planetary object as argv[2]

from numpy import *  
from numpy.fft import fft
import matplotlib.pylab as plt
import sys

# Data from files
try:
    filename1 = sys.argv[1]
    colno     = sys.argv[2]
    plan_yrs  = sys.argv[3]
except IndexError:
    filename1 = 'res2.txt'
    colno     = 5
    plan_yrs  = 12

data1 = loadtxt(filename1);

# Find no. of iterations
N1 = len(data1[:,0])

# Find max. error value
# data[x  true u(x)   est v(x)   error]
#for data in [data1, data2, data3]:
#    print("Max. error for N=%6g is:%g" % (len(data[:,0]), max(data[:,3])))
#    print("Min. error for N=%6g is:%g" % (len(data[:,0]), min(data[:,3])))


# Find radius from moon - planet: 
def radius(xcol, ycol):
    return sqrt(xcol**2 + ycol**2)
radius = vectorize(radius)

r = zeros(N1)
r = radius(data1[:,colno+2] - data1[:,colno  ],  
           data1[:,colno+3] - data1[:,colno+1])

rx = zeros(N1)
ry = zeros(N1)

rx = data1[:, colno+2] - data1[:, colno]
ry = data1[:, colno+3] - data1[:, colno+1]

rxfft =  abs(fft(rx))[:-N1/2]
ryfft =  abs(fft(ry))[:-N1/2]

fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rc({'xtick.labelsize': 20})
#plt.rc({'ytick.labelsize': 20})

ax1 = fig.add_subplot(111)
#ax2 = ax1.twinx()

#p10, = ax1.plot(data1[:,0], r, 'k-')

#p11, = ax1.plot(data1[minsindex,0]/(2*pi) , mins, 'r-')

p11, = ax1.plot(data1[:,0]/(2*pi) , rx/max(rx)) 
p12, = ax1.plot(data1[:,0]/(2*pi) , ry/max(rx))

# FFT
p13, = ax1.plot(rxfft)
p14, = ax1.plot(ryfft)

#ax1.set_xlim([10,19])
#ax1.set_xlabel('Frequency: Months/year', fontsize=14)
ax1.set_xlabel('Time [yr]', fontsize=14)
#ax1.set_ylabel(r'Rel. distance Moon-Earth [$r/r_{\rm max}$]', fontsize=14)
#ax1.set_ylim([-1, 1.4])

#plt.legend([p11,p12], \
#        [r'$x$-radius: Moon-Earth', \
#         r'$y$-radius: Moon-Earth'] )

plt.show()


