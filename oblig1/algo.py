# 2013-08-30 FYS4150
# Marius Berge Eide

# Test program implementing algorithm
#   and comparing numerical results
#   with analytical solution

from numpy import *
import matplotlib.pylab as plt

N = 100;          # no iterations
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

for i in xrange(1,N-1):
   # Forward substitution
   j = i+1
   bmm[i]  = bm[i] + (float(j-1)/j) * bmm[i-1]

bmmm[N-1]= bmm[N-1]
for i in range(N-2, 0,-1):
   # Backward substitution
   j = i+1
   bmmm[i] = bmm[i] + float(j+1)/(j+2) * bmmm[i+1]

# Normalise
for i in xrange(0,N-1):
   j = i+1
   bmmm[i] = float(j)/(j+1) * bmmm[i]


plt.figure(0)
plt.plot(x,sol, label="Analytical")
plt.plot(x,bmmm,label="Numerical")
plt.legend()
plt.show()

