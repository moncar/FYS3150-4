# 2013-08-30 FYS4150
# Marius Berge Eide

# Test program implementing algorithm
#   and comparing numerical results
#   with analytical solution

from numpy import *

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

bm = srcfunc(x)
sol= analy_func(x)

