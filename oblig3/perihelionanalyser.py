# Perihelion angle analyser module
# Reads filename holding data as argv[1]
#       first column holding x pos of object as argv[2]

from numpy import *  
import matplotlib.pylab as plt
import sys

# Data from files
try:
    filename1 = sys.argv[1]
    colno     = sys.argv[2]
    plan_yrs  = sys.argv[3]
except IndexError:
    filename1 = 'res2.txt'
    colno     = 3
    plan_yrs  = 414

data1 = loadtxt(filename1);

# Find no. of iterations
N1 = len(data1[:,0])

# Find max. error value
# data[x  true u(x)   est v(x)   error]
#for data in [data1, data2, data3]:
#    print("Max. error for N=%6g is:%g" % (len(data[:,0]), max(data[:,3])))
#    print("Min. error for N=%6g is:%g" % (len(data[:,0]), min(data[:,3])))


# Find radius from planet - Sun:
def radius(xcol, ycol):
    return sqrt(xcol**2 + ycol**2)
radius = vectorize(radius)

r = zeros(N1)
r = radius(data1[:,colno] - data1[:,1],  data1[:,colno+1] - data1[:,2])

# Find  plan_yrs+1  times perihelon angle
N2yr = N1/(plan_yrs)

mins = []
minsindex = []

min0 = min(r[0:N2yr/2])
mins.append(min0)
minsindex.append( r[0:N2yr/2].argmin())

for i in xrange(0,plan_yrs):
    currslice = r[N2yr/2 + i*N2yr : N2yr/2 + (i+1)*N2yr-1]
    mins.append(min( currslice ))
    minsindex.append(currslice.argmin() + N2yr/2 + i*N2yr)


minNN= min(r[N1-N2yr/2:-1])
mins.append(minNN)
minsindex.append( r[N1-N2yr/2:-1].argmin() + N2yr/2 + (plan_yrs-1)*N2yr)

print len(mins)

# Find perihelion angles
# tan thetaP = yP / xP

thetaP = zeros(len(minsindex))

thetaP = (data1[minsindex,colno+1]-data1[minsindex,1]) / (data1[minsindex,colno]-data1[minsindex,2])

#print thetaP



fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rc({'xtick.labelsize': 20})
#plt.rc({'ytick.labelsize': 20})

ax1 = fig.add_subplot(111)
#ax2 = ax1.twinx()

minsindex = minsindex[0:-2]
mins = mins[0:-2]

#p10, = ax1.plot(data1[:,0], r, 'k-')

print max(mins)
print min(mins)

p11, = ax1.plot(data1[minsindex,0]/(2*pi) , mins, 'r-')

plt.show()



