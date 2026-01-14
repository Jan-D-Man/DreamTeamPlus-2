import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import ScalarFormatter
from Definicions import *
from Euler import xe, ye, r_array_euler
from RK2 import x2, y2, r_array_RK2
from RK4 import x4, y4, r_array_RK4

#Error amb sol analítica

e=0.0167 #eccentricitat orbita terrestre
a=semieixMajor
c=a*e
theta_molts_euler=np.linspace(0,2*np.pi,N)
theta_molts=np.linspace(0,2*np.pi,1000)
r_molts=[]
r_molts_euler=[]

def r(theta):
  return a*(1-e**2)/(1+e*np.cos(theta))

for t in theta_molts_euler:
  r(t)
  r_molts_euler.append(r(t))

for t in theta_molts:
  r(t)
  r_molts.append(r(t))

x=np.array(r_molts)*np.cos(np.array(theta_molts))
y=np.array(r_molts)*np.sin(np.array(theta_molts))

#plt.plot(c, 0, 'ro')  # centre geomètric
#plt.plot(0, 0, 'yo') # Sol (focus)
#plt.plot(x,y)
#plt.show()

error_RK2=np.abs(np.array(r_molts)-r_array_RK2)
error_RK4=np.abs(np.array(r_molts)-r_array_RK4)
error_euler=np.abs(np.array(r_molts_euler)-r_array_euler)

# plt.plot(theta_molts,error_euler, linestyle='--')
# plt.plot(theta_molts,error_RK2)
# plt.plot(theta_molts,error_RK4)

# plt.semilogy(theta_molts, error_euler / r_molts, linestyle='--')
# plt.semilogy(theta_molts, error_RK2 / r_molts)
# plt.semilogy(theta_molts, error_RK4 / r_molts, linestyle='--')

plt.plot(theta_molts_euler, error_euler / r_molts_euler, linestyle='--')
plt.plot(theta_molts, error_RK2 / r_molts)
plt.plot(theta_molts, error_RK4 / r_molts, linestyle='--')

plt.legend(["Euler", "RK2", "RK4"])


plt.show()


