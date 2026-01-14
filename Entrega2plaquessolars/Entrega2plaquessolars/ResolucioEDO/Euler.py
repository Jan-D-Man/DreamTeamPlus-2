import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import ScalarFormatter
from Definicions import *

#euler normalitzat

theta_0=2*np.pi
DeltaTheta = 0.01/theta_0
alpha = m**2*G*M/L**2
u_0=alpha*theta_0**2
N=500
DeltaTheta = (2*np.pi/N)/theta_0


vi = 0
v_llista=[]

ui=1/(u_0*rp)
u_llista=[]
r_llista=[]
alpha1=1

for _ in range(N):

  v=vi+DeltaTheta*(alpha1-ui*(theta_0)**2)
  v_llista.append(v)

  u=ui+DeltaTheta*v

  u_llista.append(u*u_0)
  r_llista.append(1/(u_0*u))
  vi=v
  ui=u

print(v_llista,
      u_llista)

u_array = np.array(u_llista)
r_array_euler = np.array(r_llista)
theta=np.linspace(0, 2*np.pi, N)
xe=((r_array_euler)*np.cos(theta))
ye=((r_array_euler)*np.sin(theta))
phi = np.linspace(0, 2*np.pi, 1000)
xcerc= np.cos(phi)*ua
ycerc = np.sin(phi)*ua

fig, ax = plt.subplots(figsize=(8,8))
ax.grid(True,linestyle='--')
#fiquem els ticks
ax.tick_params(axis='both', which='both', direction='in', length=6)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')

#Notació científia
class CustomScalarFormatter(ScalarFormatter):
    def _set_format(self):
        self.format = "%0.2f"

formatter = CustomScalarFormatter(useMathText=True)
formatter.set_powerlimits((0, 0))
formatter.set_scientific(True)

ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_formatter(formatter)

ax.scatter(0, 0, zorder=5, label='Sol', color='yellow')

ax.plot(xe,ye, label='Euler')
ax.legend(fontsize=14)
plt.xlim(-2*10**11,2*10**11)
plt.ylim(-2*10**11,2*10**11)
#plt.plot(xcerc, ycerc, color = 'black', linestyle = '--')
plt.show()