# 2013-08-30 FYS4150
# Marius Berge Eide

# Test program implementing algorithm
#   and comparing numerical results
#   with analytical solution

from numpy import *
import matplotlib.pylab as plt

N = 300;          # no iterations
h = float(1-0)/N  # step size

x = linspace(0,1,N)

bm  = zeros(N)     # to hold b-tilde
bmm = zeros(N)     # to hold b'-tilde
bmmm= zeros(N)     # to hold b''-tilde
sol = zeros(N)     # to hold analytical solution

# Give SOURCE function f(x)
def srcfunc(x):
   return 100 * exp(-10*x)

# Give ANALYTICAL solution to Poisson eq.
def analy_func(x):
   return 1 - (1-exp(-10)) * x - exp(-10*x)
# It is tempting to vectorize and run
# both functions over x
srcfunc = vectorize(srcfunc)
analy_func = vectorize(analy_func)

bm = h**2 * srcfunc(x)
sol= analy_func(x)

# Implement algorithm
bmm[0]   = bm[0]
bmmm[N-1]= bmm[N-1]

for i in xrange(1,N-1):
   bmm[i]  = bm[i] + (float(i-1)/i) * bm[i-1]
for i in xrange(N-2, 0):
   bmmm[i] = bmm[i] + float(i+1)/(i+2) * bmm[i+1]

# Normalise
for i in xrange(0,N-1):
   bmmm[i] = float(i)/(i+1) * bmm[i]


plt.figure(0)
plt.plot(x,sol)
plt.title('Analytical')
plt.figure(1)
plt.plot(x,bmmm)
plt.title('Numerical')
plt.show()

