import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import ScalarFormatter
from Definicions import *

#paràmetres
RT, RS = 6371, 149600000
omega_T = 2 * np.pi / 24
omega_S = 2 * np.pi / (365 * 24)

# Coordenades de la UAB (Barcelona)
lat_BCN = 41.3888 * np.pi / 180
mu = (np.pi/2) - lat_BCN
delta_rad = 23.45 * np.pi / 180

# 2. M'invento els valors de les condicions inicials: t_fix és l'utilitzat a la formula de rodrigues
t_fix = np.pi / (2 * omega_T)

# Definim els dies de l'any que volem calcular
dies_per_simular = {
    'Solstici estiu (D172)': 172,
    'Solstici hivern (D355)': 355,
    'Equinocci (D80)': 80,
    'Dia 45': 45
}


fig = plt.figure(figsize=(8, 8))
ax = plt.subplot(111, polar=True)
ax.set_theta_zero_location('S') # Sud a baix
ax.set_theta_direction(-1)       # la Terra gira en sentit horari
ax.set_rlim(0, 1)


# Fem servir un mapa de colors automàtic per si volem afegir més dies
cmap = plt.get_cmap('viridis')
colors = [cmap(i) for i in np.linspace(0, 0.9, len(dies_per_simular))]

for (etiqueta, n), col in zip(dies_per_simular.items(), colors):
    # Parametrització de theta: dia 172 = pi/2
    theta_dia = (2 * np.pi * (n - 172) / 365) + (np.pi / 2)

    gamma_coords = []
    r_coords = []

    for t_rel in np.linspace(-12, 12, 2500):
        t = t_fix + t_rel

        # Vector p(t) - Rodrigues
        px = -RT * np.sin(mu) * np.sin(omega_T * t)
        py = RT * (np.cos(mu)*np.sin(delta_rad) + np.sin(mu)*np.cos(delta_rad)*np.cos(omega_T * t))
        pz = RT * (np.cos(mu)*np.cos(delta_rad) - np.sin(mu)*np.sin(delta_rad)*np.cos(omega_T * t))
        p_vec = np.array([px, py, pz])

        # Vector Sol (ara mateix està definit a partir de la parametrització (un theta diari)d'una circumferencia)
        # despres podem ferlo a partir de la solucio numerica de la edo (theta(t) per cada dia)
        rsx = -RS * np.cos(theta_dia)
        rsy = -RS * np.sin(theta_dia)
        s_vec = np.array([rsx, rsy, 0])

        # Vector distància D (rt amb l'anterior)
        d_vec = -p_vec - s_vec
        d_norm = np.linalg.norm(d_vec)

        # Angle alpha (elevació)
        sin_alpha = np.dot(d_vec, p_vec) / (d_norm * RT)
        alpha = np.arcsin(np.clip(sin_alpha, -1, 1))

        if alpha > 0: # El sol només es veu quan l'angle alpha > 0
            e_est = np.array([-np.cos(omega_T * t), -np.cos(delta_rad)*np.sin(omega_T * t), np.sin(delta_rad)*np.sin(omega_T * t)])
            e_nord = np.cross(p_vec/RT, e_est)

            d_E = np.dot(d_vec, e_est)
            d_N = np.dot(d_vec, e_nord)


            #Angle gamma (d_E per l'eix Est-Oest i -d_N per orientar el Sud a baix )
            gamma = np.arctan2(d_E, -d_N)

            gamma_coords.append(gamma)
            r_coords.append(1 - alpha/(np.pi/2))

    if gamma_coords:
        # SCATTER per evitar qualsevol línia d'unió entre punts
        ax.scatter(gamma_coords, r_coords, color=col, s=2)
        # Línia virtual per a la llegenda (sense dibuixar punts al gràfic)
        ax.plot([], [], label=etiqueta, color=col, linewidth=2)

ax.set_thetagrids([0, 90, 180, 270], labels=['Sud', 'Oest', 'Nord', 'Est'], fontsize=13)
ax.set_yticklabels([])

plt.legend(loc='upper right', frameon=True, fontsize=13)
plt.show()