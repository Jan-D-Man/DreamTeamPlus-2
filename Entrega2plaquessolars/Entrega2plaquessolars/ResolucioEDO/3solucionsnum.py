import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import ScalarFormatter
from Definicions import *
from Euler import xe, ye
from RK2 import x2, y2
from RK4 import x4, y4



#Gràfics junts

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

ax.plot(x2,y2, label='Runge-Kutta 2')
ax.plot(x4,y4, label='Runge-Kutta 4', color='blue')
ax.plot(xe,ye, label='Euler', linestyle='--',color='yellow')
ax.legend(fontsize=10, loc='upper right')
plt.xlim(-2*10**11,2*10**11)
plt.ylim(-2*10**11,2*10**11)
#plt.plot(xcerc, ycerc, color = 'black', linestyle = '--')
plt.show()