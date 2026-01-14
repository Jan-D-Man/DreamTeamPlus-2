import numpy as np
import matplotlib.pyplot as plt
from Definicions import *
from ResolucioEDO.RK4 import r_array_RK4, DeltaTheta


# =============================================================================
# 1. PARÀMETRES BASE (UAB)
# =============================================================================
RT, RS = 6371, 149600000
omega_T = 2 * np.pi / 24
lat_BCN = 41.3888 * np.pi / 180
mu = (np.pi/2) - lat_BCN
delta_rad = 23.45 * np.pi / 180
I_0, area, eficiencia, p_max = 1381, 2, 0.2, 400

dies_any = np.arange(1, 366)
# Augmentem resolució temporal dins del dia per a una integral més precisa
hores_dia = np.linspace(0, 24, 200)
cmap = plt.get_cmap('tab10')

def simulacio_vectoritzada(b_graus, g_graus):
    b_rad = np.radians(b_graus)
    g_rad = np.radians(g_graus)
    mitjanes_diaries = []

    for n in dies_any:
        theta_dia = (2 * np.pi * (n - 172) / 365) - (np.pi / 2)
        t_central = (theta_dia / omega_T) + (np.pi / (2 * omega_T))

        # Vectorització de les hores del dia
        t = t_central + (hores_dia - 12)
        arg = omega_T * t

        p_unit = np.array([
            -np.sin(mu) * np.sin(arg),
            (np.cos(mu)*np.sin(delta_rad) + np.sin(mu)*np.cos(delta_rad)*np.cos(arg)),
            (np.cos(mu)*np.cos(delta_rad) - np.sin(mu)*np.sin(delta_rad)*np.cos(arg))
        ])

        s_unit = np.array([np.cos(theta_dia), np.sin(theta_dia), 0])
        d_unit = -s_unit # El Sol està tan lluny que els rajos són paral·lels

        # Elevació alpha i màscara de dia
        alpha = np.arcsin(np.clip(np.dot(d_unit, p_unit), -1, 1))
        mask_sol = alpha > 0

        potencies = np.zeros_like(hores_dia)

        if np.any(mask_sol):
            p_sol = p_unit[:, mask_sol]
            arg_sol = arg[mask_sol]

            # Base local
            e_est = np.array([-np.cos(arg_sol), -np.cos(delta_rad)*np.sin(arg_sol), np.sin(delta_rad)*np.sin(arg_sol)])
            e_nord = np.cross(p_sol.T, e_est.T).T

            # Normal de la placa
            n_placa = (np.cos(b_rad) * p_sol +
                       np.sin(b_rad) * (-np.cos(g_rad) * e_nord + np.sin(g_rad) * e_est))

            cos_inc = np.einsum('i,ij->j', d_unit, n_placa)
            pot_sol = area * eficiencia * I_0 * np.maximum(0, cos_inc)
            pot_sol = np.minimum(p_max, pot_sol)
            potencies[mask_sol] = pot_sol

        # Energia del dia = Integral de la potència (Wh)
        energia_dia_Wh = np.trapz(potencies, x=hores_dia)
        # Potència mitjana = Energia / 24 hores
        mitjanes_diaries.append(energia_dia_Wh / 24)

    return np.array(mitjanes_diaries)

# =============================================================================
# 3. GENERACIÓ DE GRÀFIQUES COMPARATIVES
# =============================================================================
plt.figure(figsize=(15, 9))

# --- GRÀFICA 1: VARIACIÓ D'ORIENTACIÓ ---
plt.subplot(2, 1, 1)
for i, g in enumerate([-90, -30, 0, 30, 90]):
    mitjanes = simulacio_vectoritzada(45, g)
    # L'energia total és la integral de la potència mitjana al llarg del temps (en hores)
    # Com que l'eix X són dies (1 unitat = 24 hores), multipliquem per 24
    energia_anual_kWh = np.trapz(mitjanes, dx=1) * 24 / 1000
    plt.plot(dies_any, mitjanes, color=cmap(i), label=f"Azimut {g} (E={energia_anual_kWh:.1f} kWh)")

plt.ylabel("Potencia Mitjana (W)", fontsize=12)
plt.xlabel("Dia de l'any", fontsize=12)
plt.legend(loc='upper right', fontsize=13)
plt.xlim(1,400)
plt.ylim(20,200)
plt.grid(alpha=0.3)
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