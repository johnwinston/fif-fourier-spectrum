import cmath
import scipy
import quadpy
import numpy as np
import matplotlib.pyplot as plt

values = [0, 0, 0.25, 1, 0.5, 1.4, 0.75, -0.5, 1, 0]
timeSteps = 1000
globalIt = 4

d = [-.5, .7, .7, .7]
a = []
e = []
c = []
f = []

nMaps = int((len(values)/2)-1)
print("nMaps: ", nMaps)

# FIF
def FIF(x, it):
   if (it == 0):
      p = np.array([(c[i]*((x - e[i]) / a[i]) + f[i]) * uFunc(x,i)
                     for i in range(nMaps)])
      return p.sum()
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

###########################################################

time = np.linspace(0, 1,num=timeSteps)
simpleParameters(values)

FIFPoints = np.array([FIF(t, globalIt) for t in time])
plt.plot(time, FIFPoints)

print(values)
print("a:", a)
print("e:", e)
print("c:", c)
print("f:", f)
print("d:", d)
plt.show()
