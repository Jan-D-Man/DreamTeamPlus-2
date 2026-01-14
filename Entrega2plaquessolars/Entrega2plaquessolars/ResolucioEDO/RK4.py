import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import ScalarFormatter
from Definicions import *


#RK 4 NORMALITZAT
theta_0=2*np.pi
alpha = m**2*G*M/L**2
u_0=alpha*theta_0**2
N=1000
DeltaTheta = (2*np.pi/N)/theta_0

ui=1/(u_0*rp)
vi=0
alpha1=1
u_llista=[]
r_llista=[]

for _ in range(N):

    # k1
    k1_u=vi
    k1_v=alpha1-ui*theta_0**2

    # k2
    k2_u=vi+0.5*DeltaTheta*k1_v
    k2_v=alpha1-(ui+0.5*DeltaTheta*k1_u)*theta_0**2

    # k3
    k3_u=vi+0.5*DeltaTheta*k2_v
    k3_v=alpha1-(ui+0.5*DeltaTheta*k2_u)*theta_0**2

    # k4
    k4_u=vi+DeltaTheta*k3_v
    k4_v=alpha1-(ui+DeltaTheta*k3_u)*theta_0**2

    # actualització
    ui =ui+(DeltaTheta/6)*(k1_u+2*k2_u+2*k3_u+k4_u)
    vi =vi+(DeltaTheta/6)*(k1_v+2*k2_v+2*k3_v+k4_v)

    u_llista.append(ui*u_0)
    r_llista.append(1/(u_0*ui))

u_array = np.array(u_llista)
r_array_RK4 = np.array(r_llista)
theta=np.linspace(0, 2*np.pi, N)
x4=((r_array_RK4)*np.cos(theta))
y4=((r_array_RK4)*np.sin(theta))
phi = np.linspace(0, 2*np.pi, 1000)
xcerc= np.cos(phi)*ua
ycerc = np.sin(phi)*ua

fig, ax2 = plt.subplots(figsize=(8,8))
ax2.grid(True,linestyle='--')
#fiquem els ticks
ax2.tick_params(axis='both', which='both', direction='in', length=6)
ax2.xaxis.set_ticks_position('both')
ax2.yaxis.set_ticks_position('both')

#Notació científia
class CustomScalarFormatter(ScalarFormatter):
    def _set_format(self):
        self.format = "%0.2f"

formatter = CustomScalarFormatter(useMathText=True)
formatter.set_powerlimits((0, 0))
formatter.set_scientific(True)

ax2.xaxis.set_major_formatter(formatter)
ax2.yaxis.set_major_formatter(formatter)

ax2.scatter(0, 0, zorder=5, label='Sol', color='yellow')

ax2.plot(x4,y4, label='Runge-Kutta 4')
ax2.legend(fontsize=14)
plt.xlim(-2*10**11,2*10**11)
plt.ylim(-2*10**11,2*10**11)
#plt.plot(xcerc, ycerc, color = 'black', linestyle = '--')
plt.show()