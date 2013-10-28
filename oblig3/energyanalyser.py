# Kinetic energy analyser module
# Reads filename holding data as argv[1]
#       first column holding x pos of planetary object as argv[2]

from numpy import *  
from numpy.fft import fft
import matplotlib.pylab as plt
import sys

# Constants
fMsol = 332948.6;   # One Solar mass in Earth masses

# Data from files
try:
    filename_pos = sys.argv[1]
    filename_vel = sys.argv[2]
    filename_energy=sys.argv[3]
    colno     = sys.argv[4]
    mass      = sys.argv[5]
except IndexError:
    filename_pos = 'res2.txt'
    filename_vel = 'res2vel.txt'
    filename_energy= 'res2energy.txt'
    colno     = 3
    mass      = 0.815/fMsol # Mercury
    print "Applying default settings for celestial object: column ",colno


#data1 = loadtxt(filename_pos);
data2 = loadtxt(filename_vel);
data3 = loadtxt(filename_energy);

# Find no. of iterations
N1 = len(data2[:,0])

# Find max. error value
# data[x  true u(x)   est v(x)   error]
#for data in [data1, data2, data3]:
#    print("Max. error for N=%6g is:%g" % (len(data[:,0]), max(data[:,3])))
#    print("Min. error for N=%6g is:%g" % (len(data[:,0]), min(data[:,3])))


# Find kinetic energy of celestial object located at x: colno, y: colno+1
def kinetic(vx,vy,mass):
    return 0.5 * mass * (vx*vx + vy*vy)
kinetic = vectorize(kinetic)

kinenergy = zeros(N1)
kinenergy = kinetic(vx=data2[:,colno], vy=data2[:,colno+1], mass=mass)

totVirial = 2*kinenergy + data3[:,1]

meanval = mean(abs(data3[:,1])) 

# Scale energy: to initial value, which is the value it should retain

fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rc({'xtick.labelsize': 20})
#plt.rc({'ytick.labelsize': 20})

ax1 = fig.add_subplot(111)
#ax2 = ax1.twinx()


p10, = ax1.plot(data2[:,0]/(2*pi) , kinenergy/meanval)
p11, = ax1.plot(data3[:,0]/(2*pi) , data3[:,1]/meanval)
p12, = ax1.plot(data3[:,0]/(2*pi) , totVirial/meanval)


ax1.set_xlabel('Time [yr]', fontsize=14)
ax1.set_ylabel(r'Rel. energy [$E/ \langle |V| \rangle $]', fontsize=14)
ax1.set_ylim([-1.2, 0.6])

plt.legend([p10,p11,p12], \
        [r'Kinetic energy $K$', \
         r'Potential energy $V$', 
         r'Virial theorem: $2\langle K \rangle + \langle V \rangle = 0$'],
        loc=5)

plt.show()

