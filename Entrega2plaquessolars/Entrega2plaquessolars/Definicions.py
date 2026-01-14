import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import ScalarFormatter


ua = 149597870700
rp= 147098290*10**3
ra = 152098232*10**3
semieixMajor = ua
semieixMenor = 0.99986*ua
G = 6.67428*10**-11
M = 1.98847*10**30
m = 5.97*10**24
s = np.pi*semieixMajor*semieixMenor
T = 365*24*3600
L = 2*m*s/T

DeltaTheta = 0.05
alpha = m**2*G*M/L**2
N = int(2 * np.pi / DeltaTheta)