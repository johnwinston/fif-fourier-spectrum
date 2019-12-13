import cmath
import scipy
import quadpy
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from matplotlib.pyplot import draw, show

values = [0, 0, 0.25, 1, 0.5, 1.4, 0.75, -0.5, 1, 0]
timeSteps = 1000
globalIt = 0

d = [-.5, .8, .8, .8]
a = []
e = []
c = []
f = []

nMaps = int((len(values)/2)-1)
print("nMaps: ", nMaps)

# Picecewise linear function between interpolation points
def pieceWise(x):
   z = np.array([(c[i]*((x - e[i]) / a[i]) + f[i]) * uFunc(x,i)
                 for i in range(nMaps)])
   return z.sum()

# FIF
def FIF(x, it):
   if (it == 0):
      l = np.array([(d[i]*(pieceWise((x - e[i]) / a[i])) +
                     c[i]*((x - e[i]) / a[i]) + f[i]) * uFunc(x,i)
                     for i in range(nMaps)])
      return l.sum()
   p = np.array([(d[i]*(FIF((x - e[i]) / a[i], it-1)) +
                  c[i]*((x - e[i]) / a[i]) + f[i]) * uFunc(x,i)
                  for i in range(nMaps)])
   return p.sum()

# Returns 1 or 0 based on section between interpolation points x is in
def uFunc(x,i):
   if (i+1 > nMaps-1):
      return 1 if x >= e[i] else 0
   u1 = 1 if x - e[i] >= 0 else 0
   u2 = 1 if x - e[i+1] >= 0 else 0
   return u1 - u2

# Used when FIF is expected to be between 0 and 1 with evenly space points
def simpleParameters(values):
   k = 2
   for i in range(nMaps):
      a.append(1 / nMaps)
      e.append(i / nMaps)
      c.append(values[k+1] - values[k-1])
      f.append(values[k-1])
      k=k+2

def gOmega(x):
   if (x == 0):
      return 0
   p = np.array([(cmath.exp((-1j * x) * (i / nMaps))  -
                  cmath.exp((-1j * x) * ((i + 1) / nMaps)))
                  / ((1j * x) * (1j * x)) 
                 for i in range(nMaps)])
   o = np.array([((c[i] / a[i]) + ((1j * x) *
                 (f[i] - ((c[i]*e[i])/a[i]))))
                 for i in range(nMaps)])
   l = np.array([(c[i] / a[i]) * ((((i / nMaps) *
                 (cmath.exp((-1j * x) * (i / nMaps)))) -
                 (((i + 1) / nMaps) * (cmath.exp((-1j * x) *
                 ((i + 1) / nMaps))))) / (1j * x))
                 for i in range(nMaps)])
   k = np.array([(p[i] * o[i]) + l[i]
                 for i in range(nMaps)])
   return np.sqrt(np.power(scipy.real(k.sum()), 2) +
                  np.power(scipy.imag(k.sum()),2))

# DFT of vertical scaling factors
def qOmega(x):
   y = np.array([d[i]*cmath.exp(-1j*((x*i)/nMaps))
                 for i in range(nMaps)])
   return np.sqrt(np.power(scipy.real(a[0]*y.sum()), 2) +
                  np.power(scipy.imag(a[0]*y.sum()),2))

# Spectrum
def quadInt(func):
   def real_func(x):
      f = FIF(x, globalIt)
      return scipy.real(f*func(x))
   def imag_func(x):
      f = FIF(x, globalIt)
      return scipy.imag(f*func(x))
   real_integral = quad(real_func, 0, 1)
   imag_integral = quad(imag_func, 0, 1)
   return np.sqrt(np.power(real_integral[0], 2) +
                  np.power(imag_integral[0],2))

def everything(t):
   return qOmega(t)*quadInt(lambda x: (scipy.exp(-1j*(a[0]*t)*x)))+gOmega(t)

###########################################################

time = np.linspace(0, 1,num=timeSteps)
simpleParameters(values)

FIFPoints = np.array([FIF(t, globalIt) for t in time])
plt.subplot(211)
plt.plot(time, FIFPoints)

time = np.linspace(1, 100, num=timeSteps)
integrationPoints = np.array([everything(t) for t in time])
plt.subplot(212)
plt.plot(time, integrationPoints)

print(values)
print("a:", a)
print("e:", e)
print("c:", c)
print("f:", f)
print("d:", d)
plt.show()
