# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from Definicions import *
from ResolucioEDO.RK4 import r_array_RK4, DeltaTheta



RT, RS = 6371, 149600000
omega_T = 2 * np.pi / 24
lat_BCN = 41.3888 * np.pi / 180
mu = (np.pi/2) - lat_BCN
delta_rad = 23.45 * np.pi / 180
# Calculem el temps total de l'òrbita
temps = int(sum(r_array_RK4**2*DeltaTheta)*m/L)
TEMPS = []
z = 0
# calculam els temps corresponents radi amb el sumatori acumulat per a cada cas
for r_i in r_array_RK4:
    t_i = sum(r_array_RK4[:z]**2*DeltaTheta)*m/L
    TEMPS.append(t_i)
    z = z +1

# Transformem el tems a dies
dies = np.array(TEMPS)/(3600*24)
dies = np.round(np.array(TEMPS)/(3600*24))
r_dies = np.array([r_array_RK4[np.round(TEMPS) == d].mean() for d in dies])
r_dies = np.roll(r_dies, -172)

I_0, beta, gamma_p, area, eficiencia_total, p_max = 1381, 0, 0, 2, 0.2, 400


dies_totals = np.arange(1, 366)
hores_per_dia = np.linspace(0, 24, 100) # Més resolució per a una integral precisa

temps_eix_x = []
potencia_eix_y = []
mitjanes_diaries = []
dies_eix_mitjanes = []

for n in dies_totals:
    theta_dia = (2 * np.pi * (n - 172) / 365) - (np.pi / 2)
    t_central_fisic = (theta_dia / omega_T) + (np.pi / (2 * omega_T))

    potencies_avui = []

    for hora_lst in hores_per_dia:
        t_fisic = t_central_fisic + (hora_lst - 12)
        arg = omega_T * t_fisic

        # Vectors
        px = -RT * np.sin(mu) * np.sin(arg)
        py = RT * (np.cos(mu)*np.sin(delta_rad) + np.sin(mu)*np.cos(delta_rad)*np.cos(arg))
        pz = RT * (np.cos(mu)*np.cos(delta_rad) - np.sin(mu)*np.sin(delta_rad)*np.cos(arg))
        p_vec = np.array([px, py, pz])
        p_unit = p_vec / RT

        s_vec = np.array([r_dies[n] * np.cos(theta_dia), r_dies[n] * np.sin(theta_dia), 0])
        d_vec = -p_vec - s_vec
        d_norm = np.linalg.norm(d_vec)
        d_unit = d_vec / d_norm

        sin_alpha = np.dot(d_vec, p_vec) / (d_norm * RT)    
        alpha = np.arcsin(np.clip(sin_alpha, -1, 1))

        pot = 0
        if alpha > 0:
            e_est = np.array([-np.cos(arg), -np.cos(delta_rad)*np.sin(arg), np.sin(delta_rad)*np.sin(arg)])
            e_nord = np.cross(p_unit, e_est)
            n_placa = (np.cos(beta) * p_unit +
                       np.sin(beta) * (-np.cos(gamma_p) * e_nord + np.sin(gamma_p) * e_est))

            cos_incidencia = np.dot(d_unit, n_placa)
            if cos_incidencia > 0:
                pot = area * eficiencia_total * I_0 * cos_incidencia
                if pot > p_max: pot = p_max

        potencies_avui.append(pot)
        temps_eix_x.append(n + (hora_lst / 24))
        potencia_eix_y.append(pot)

    # Càlcul de la mitjana d'avui (Integral del dia / 24)
    energia_avui_Wh = np.trapz(potencies_avui, x=hores_per_dia)
    mitjanes_diaries.append(energia_avui_Wh / 24)
    dies_eix_mitjanes.append(n + 0.5) # Centrem el punt al mig del dia

# Energia total anual integrant les mitjanes (o sumant l'energia diària)
energia_total_kWh = np.sum(mitjanes_diaries) * 24 / 1000


fig, ax = plt.subplots(figsize=(15, 5))

# 1. Les 365 corbes (Potència instantània)
ax.plot(temps_eix_x, potencia_eix_y, color='teal', linewidth=0.3, alpha=0.4, label='Potencia instantania')

# 2. La línia de potència mitjana diària
ax.plot(dies_eix_mitjanes, mitjanes_diaries, color='red', linewidth=2, label='Potencia mitjana diaria')

# Estètica
mesos_noms = ['Gen', 'Feb', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Oct', 'Nov', 'Des']
mesos_dies = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
ax.set_xticks(mesos_dies)
ax.set_xticklabels(mesos_noms)

ax.set_xlim(1, 365)
ax.set_ylim(0, 450)
ax.set_xlabel("Mes", fontsize=13)
ax.set_ylabel("Potencia (W)", fontsize=13)
ax.legend(loc='upper right', title_fontsize=13, fontsize=13)
ax.grid(True, linestyle=':', alpha=0.6)


plt.tick_params(
    axis='both',
    direction='in',       # Guionets cap a dins
    top=True,             # Guionets al costat superior
    right=True,           # Guionets al costat dret
    which='both',
    labelsize=12
)


plt.tight_layout()
plt.show()