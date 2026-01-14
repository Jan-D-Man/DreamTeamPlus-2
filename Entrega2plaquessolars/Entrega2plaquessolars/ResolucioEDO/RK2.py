import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import ScalarFormatter
from Definicions import *

#RK 2 NORMALITZAT

theta_0=2*np.pi
alpha = m**2*G*M/L**2
u_0=alpha*theta_0**2
N=1000
DeltaTheta = (2*np.pi/N)/theta_0


vi = 0
v_llista=[]

ui=1/(u_0*rp)
u_llista=[]
r_llista=[]
alpha1=1

ui=1/(rp*u_0)
vi=0

u_llista = []
r_llista = []

for _ in range(N):

    # k1
    k1_u = vi
    k1_v = alpha1 - ui*theta_0**2

    # punt mig
    u_mid=ui+0.5*DeltaTheta*k1_u
    v_mid=vi+0.5*DeltaTheta*k1_v

    # k2
    k2_u=v_mid
    k2_v=alpha1-u_mid*theta_0**2

    # actualització
    ui= ui+DeltaTheta*k2_u
    vi= vi+DeltaTheta*k2_v

    u_llista.append(ui*u_0)
    r_llista.append(1/(ui*u_0))

u_array = np.array(u_llista)
r_array_RK2 = np.array(r_llista)
theta=np.linspace(0, 2*np.pi, N)
x=((r_array_RK2)*np.cos(theta))
y=((r_array_RK2)*np.sin(theta))
x2=((r_array_RK2)*np.cos(theta))
y2=((r_array_RK2)*np.sin(theta))
phi = np.linspace(0, 2*np.pi, 1000)
xcerc= np.cos(phi)*ua
ycerc = np.sin(phi)*ua
#print(np.cos(theta))
#print(x, "sexe", y)



fig, ax1 = plt.subplots(figsize=(8,8))
ax1.grid(True,linestyle='--')
#fiquem els ticks
ax1.tick_params(axis='both', which='both', direction='in', length=6)
ax1.xaxis.set_ticks_position('both')
ax1.yaxis.set_ticks_position('both')

#Notació científia
class CustomScalarFormatter(ScalarFormatter):
    def _set_format(self):
        self.format = "%0.2f"

formatter = CustomScalarFormatter(useMathText=True)
formatter.set_powerlimits((0, 0))
formatter.set_scientific(True)

ax1.xaxis.set_major_formatter(formatter)
ax1.yaxis.set_major_formatter(formatter)

ax1.scatter(0, 0, zorder=5, label='Sol', color='yellow')

ax1.plot(x2,y2, label='Runge-Kutta 2')
ax1.legend(fontsize=14)
plt.xlim(-2*10**11,2*10**11)
plt.ylim(-2*10**11,2*10**11)
#plt.plot(xcerc, ycerc, color = 'black', linestyle = '--')
plt.show()


