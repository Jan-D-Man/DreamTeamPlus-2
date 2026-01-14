import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import ScalarFormatter

# Paràmetres globals 
ua = 149597870700 # unitat astronòmica en metres
rp= 147098290*10**3 # radi periheli en metres
ra = 152098232*10**3 # radi afeli en metres
semieixMajor = ua # semieix major en metres
semieixMenor = 0.99986*ua # semieix menor en metres
G = 6.67428*10**-11 # constant gravitacional en m^3 kg^-1 s^-2
M = 1.98847*10**30 # massa del Sol en kg
m = 5.97*10**24 # massa de la Terra en kg
s = np.pi*semieixMajor*semieixMenor # àrea de l'òrbita terrestre en m^2
T = 365*24*3600 # període orbital de la Terra en segons
L = 2*m*s/T # moment angular específic de la Terra en m^2 kg s^-1
RT, RS = 6371, 149600000 # radi de la Terra i distància mitjana al Sol en km
DeltaTheta = 0.001 # pas en radians per a la resolució de l'EDO
alpha = m**2*G*M/L**2 # paràmetre de l'òrbita
N = int(2 * np.pi / DeltaTheta) # nombre de passos per a una òrbita completa

omega_T = 2 * np.pi / 24 # velocitat angular de rotació de la Terra en rad/h
lat_BCN = 41.3888 * np.pi / 180 # latitud de Barcelona en radians
mu = (np.pi/2) - lat_BCN # complement de la latitud
delta_rad = 23.45 * np.pi / 180 # obliqüitat de l'eclíptica en radians

I_0, area, eficiencia, p_max = 1381, 2, 0.2, 400 # irradiància solar (W/m^2), àrea placa (m^2), eficiència i potència màxima (W)