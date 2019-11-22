import numpy as np
import scipy, random
import matplotlib.pyplot as plt
from matplotlib.pyplot import draw, show
import math
import cmath

from scipy.integrate import quad
import quadpy

fig, ax = plt.subplots()

values = [0, 0, 0.25, 1, 0.5, 1.4, 0.75, -0.5, 1, 0]
gOmegaPoints = []
qOmegaPoints = []
integrationPoints = []
integrationPointsFull = []
FIFPoints = []

globalIt = 4

d = .4
a = []
e = []
c = []
f = []

nMaps = int((len(values)/2)-1)
sum = 0

def calcValues(values):
   del a[:]
   del e[:]
   del c[:]
   del f[:]
   k = 2
   for i in range((len(values)/2)-1):
      a.append((values[k] - values[k-2]) / float(values[len(values)-2] - values[0]))
      e.append((values[len(values)-2]*values[k-2] - values[0]*values[k]) / 
         float(values[len(values)-2] - values[0]))
      c.append(((values[k+1] - values[k-1]) / float(values[len(values)-2] - values[0]))
         - ((d*(values[len(values)-1] - values[1])) / float(values[len(values)-2] - values[0])))
      f.append(((values[len(values)-2]*values[k-1] - values[0]*values[k+1])
         / float(values[len(values)-2]-values[0])) - \
      ((d*(values[len(values)-2]*values[1] - values[0]*values[len(values)-1]))
         / float(values[len(values)-2] - values[0])))
      k = k + 2

def plotPoints(values):
   k = 0
   for i in range(int(len(values) / 2)):
      plt.plot(values[k],values[k+1], '.',markersize=2)
      k = k + 2
	
# Picecewise linear function between interpolation points
def pieceWise(x):
   sum = 0
   for i in range(nMaps):
      sum += (c[i]*((x - e[i]) / float(a[i])) + f[i]) * uFunc(x,i)
   return sum

# Non-recursive FIF
def simpleFIF(x):
   sum = 0
   for i in range(nMaps):
      sum += (d*(pieceWise((x - e[i]) / float(a[i]),nMaps)) +
         c[i]*((x - e[i]) / float(a[i])) + f[i]) * uFunc(x,i)
   gOmegaPoints.append(x)
   gOmegaPoints.append(sum)
   return sum
   
# Recursive FIF
def complexFIF(x, it):
   sum = 0
   if (it == 0):
      for i in range(nMaps):
         sum += (d*(pieceWise((x - e[i]) / float(a[i])))
            + c[i]*((x - e[i]) / float(a[i])) + f[i]) * uFunc(x,i)
      return sum
	
   for i in range(nMaps):
      sum += (d*(complexFIF((x - e[i]) / float(a[i]), it-1))
         + c[i]*((x - e[i]) / float(a[i])) + f[i]) * uFunc(x,i)

   return sum

# Returns 1 or 0 based on section between interpolation points x is in
def uFunc(x,i):

   if (i+1 > nMaps-1):
      if(x >= e[i]):
         return 1
      elif(x < e[i]):
         return 0
   if((x - e[i]) >= 0):
      u1 = 1
   elif((x - e[i]) < 0):
      u1 = 0
	
   if((x - e[i+1]) >= 0):
      u2 = 1
   elif((x - e[i+1]) < 0):
      u2 = 0
   return float(u1 - u2)

# Used when FIF is expected to be between 0 and 1 with evenly space points
def simpleParameters(values):
   del a[:]
   del e[:]
   del c[:]
   del f[:]
   k = 2
   for i in range((len(values)/2)-1):
      a.append(1 / float(nMaps))
      e.append(i / float(nMaps))
      c.append(values[k+1] - values[k-1])
      f.append(values[k-1])
      k=k+2
	  

def gOmega(x):
   gsum = 0
   for i in range(nMaps - 1):
      first = (cmath.exp(-1j * x * (i / float(nMaps)))  - 
         cmath.exp(-1j * x * ((i + 1) / float(nMaps))))  / ((1j * x) * (1j * x))
      second = ((c[i] / float(a[i])) + ((1j * x) * (f[i] - ((c[i]*e[i])/float(a[i])))))
 
      third = (c[i] / float(a[i])) * (((((i / float(nMaps)) * (cmath.exp(-1j * x * (i / float(nMaps))))))
         - (((i + 1) / float(nMaps)) * (cmath.exp(-1j * x * ((i + 1) / float(nMaps)))))) / (1j * x))
      gsum += ((first * second) + third)

   return np.sqrt(np.power(scipy.real(gsum), 2) + np.power(scipy.imag(gsum),2))


def qOmega(x):
   sum = 0
   for i in range(nMaps - 1):
      sum += d*cmath.exp(-1j*x*(i/float(nMaps)))
   sum = sum * a[0]
   return np.sqrt(np.power(scipy.real(sum), 2) + np.power(scipy.imag(sum),2))
   
def quadInt(func):

   def real_func(x):
      f = complexFIF(x, globalIt, 0)
      return scipy.real(f*func(x))
   
   def imag_func(x):
      f = complexFIF(x, globalIt, 0)
      return scipy.imag(f*func(x))

   real_integral = quad(real_func, 0, 1)
   imag_integral = quad(imag_func, 0, 1)
   return np.sqrt(np.power(real_integral[0], 2) + np.power(imag_integral[0],2))

def everything(t):
   return np.sqrt(np.power(scipy.real(qOmega(t)*quadInt(lambda x: (scipy.exp(-1j*(a[0]*t)*x)))+gOmega(t)), 2) + np.power(scipy.imag(qOmega(t)*quadInt(lambda x: (scipy.exp(-1j*(a[0]*t)*x)))+gOmega(t)), 2))

########################################################### Begin script
#simpleParameters(values)
# For interactive graph
#plt.ion()
#plt.show(block=False)
a = [0.25, 0.25, 0.25, 0.25]
e = [0.0, 0.25, 0.5, 0.75]
c = [0.7, .075, -.93, .155]
f = [0, .7, .775, -.155]

'''
hz = 0
rangeI = 40
sum = 0
for i in range(rangeI):
   sum += (1 / float(rangeI))
   hz += 1/float(2)
   qOmegaPoints.append(hz)
   qOmegaPoints.append(qOmega(hz))
   gOmegaPoints.append(hz)
   gOmegaPoints.append(gOmega(hz))
   integrationPoints.append(hz)
   integrationPoints.append(quadInt(lambda x: (scipy.exp(-1j*hz*(a[0]*x)))))
   print "Percent Complete:", sum * 100
'''

hz = 0
rangeI = 1000
sum = 0

# FIF looks broken
for i in range(rangeI):
   sum += (1 / float(rangeI))
   hz += 1/float(rangeI)
   #integrationPointsFull.append(hz)
   #integrationPointsFull.append(everything(hz))
   #integrationPoints.append(hz)
   #integrationPoints.append(quadInt(lambda x: (scipy.exp(-1j*hz*x))))
   FIFPoints.append(hz)
   FIFPoints.append(pieceWise(hz))
   integrationPoints.append(hz)
   integrationPoints.append(complexFIF(hz, globalIt))
   print("Percent Complete:", sum * 100)

plt.subplot(211)
plotPoints(FIFPoints)
plt.subplot(212)
plotPoints(integrationPoints)
'''
plt.subplot(4,1,1)
plotPoints(qOmegaPoints)
plt.subplot(4,1,2)
plotPoints(gOmegaPoints)
plt.subplot(4,1,3)
plotPoints(integrationPoints)
plt.subplot(4,1,4)
plotPoints(integrationPointsFull)
'''
print(values)
print("a:", a)
print("e:", e)
print("c:", c)
print("f:", f)
print("d:", d)
plt.show(block=True)
