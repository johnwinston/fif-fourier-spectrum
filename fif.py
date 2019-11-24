import cmath
import scipy
import quadpy
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from matplotlib.pyplot import draw, show

values = [0, 0, 0.25, 1, 0.5, 1.4, 0.75, -0.5, 1, 0]
timeSteps = 1000
globalIt = 2

d = [.4, .8, .7, .4]
a = []
e = []
c = []
f = []

nMaps = (len(values)/2)-1
	
# Picecewise linear function between interpolation points
def pieceWise(x):
   sum = 0
   for i in range(int(nMaps)):
      sum += (c[i]*((x - e[i]) / float(a[i])) + f[i]) * uFunc(x,i)
   return sum
   
# FIF
def FIF(x, it):
   sum = 0
   if (it == 0):
      for i in range(int(nMaps)):
         sum += (d[i]*(pieceWise((x - e[i]) / float(a[i])))
            + c[i]*((x - e[i]) / float(a[i])) + f[i]) * uFunc(x,i)
      return sum
	
   for i in range(int(nMaps)):
      sum += (d[i]*(FIF((x - e[i]) / float(a[i]), it-1))
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
   for i in range((int(len(values)/2)-1)):
      a.append(1 / float(nMaps))
      e.append(i / float(nMaps))
      c.append(values[k+1] - values[k-1])
      f.append(values[k-1])
      k=k+2

def gOmega(x):
   gsum = 0
   for i in range(int(nMaps)):
      first = (cmath.exp((-1j * x) * (i / float(nMaps)))  - 
         cmath.exp((-1j * x) * ((i + 1) / float(nMaps)))) / ((1j * x) * (1j * x))
      second = ((c[i] / float(a[i])) + ((1j * x) * (f[i] - ((c[i]*e[i])/float(a[i])))))
 
      third = (c[i] / float(a[i])) * ((((i / float(nMaps)) * (cmath.exp((-1j * x) * (i / float(nMaps)))))
         - (((i + 1) / float(nMaps)) * (cmath.exp((-1j * x) * ((i + 1) / float(nMaps)))))) / (1j * x))
      gsum += ((first * second) + third)
   return np.sqrt(np.power(scipy.real(gsum), 2) + np.power(scipy.imag(gsum),2))

def qOmega(x):
   sum = 0
   for i in range(int(nMaps)):
      sum += d[i]*cmath.exp(-1j*((x*i)/float(nMaps)))
   sum = sum * a[0]
   return np.sqrt(np.power(scipy.real(sum), 2) + np.power(scipy.imag(sum),2))
   
def quadInt(func):
   def real_func(x):
      f = FIF(x, globalIt)
      return scipy.real(f*func(x))
   def imag_func(x):
      f = FIF(x, globalIt)
      return scipy.imag(f*func(x))
   real_integral = quad(real_func, 0, 1)
   imag_integral = quad(imag_func, 0, 1)
   return np.sqrt(np.power(real_integral[0], 2) + np.power(imag_integral[0],2))

def everything(t):
   return qOmega(t)*quadInt(lambda x: (scipy.exp(-1j*(a[0]*t)*x)))+gOmega(t)
   
def iterative(x, its):
   sum = 0
   sumTemp = 0
   mSum = 0
   for i in range(its):
      for j in range(i-1):
        if (sumTemp != 0):
           sumTemp *= qOmega(x*(np.power(a[0],j)))
        else:
           sumTemp = qOmega(x*(np.power(a[0],j)))
      sum += gOmega(x*(np.power(a[0], i))) * sumTemp
      sumTemp = 0
   return sum

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
