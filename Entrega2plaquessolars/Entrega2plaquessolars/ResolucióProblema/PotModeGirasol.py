import numpy as np
import matplotlib.pyplot as plt
from Definicions import *


dies_any = np.arange(1, 366)
hores_dia = np.linspace(0, 24, 250)

def simulacio_seguidors(tipus="fix", b_fix=36):
    mitjanes_diaries = []

    for n in dies_any:
        theta_dia = (2 * np.pi * (n - 172) / 365) - (np.pi / 2)
        t_central = (theta_dia / omega_T) + (np.pi / (2 * omega_T))

        t = t_central + (hores_dia - 12)
        arg = omega_T * t

        p_unit = np.array([
            -np.sin(mu) * np.sin(arg),
            (np.cos(mu)*np.sin(delta_rad) + np.sin(mu)*np.cos(delta_rad)*np.cos(arg)),
            (np.cos(mu)*np.cos(delta_rad) - np.sin(mu)*np.sin(delta_rad)*np.cos(arg))
        ])

        s_unit = np.array([np.cos(theta_dia), np.sin(theta_dia), 0])
        d_unit = -s_unit

        # Elevació alpha per saber si és de dia
        cos_alpha = np.einsum('i,ij->j', d_unit, p_unit)
        alpha = np.arcsin(np.clip(cos_alpha, -1, 1))
        mask_sol = alpha > 0

        potencies = np.zeros_like(hores_dia)

        if np.any(mask_sol):
            if tipus == "2-eixos":
                # Cas 1: Cos_inc sempre és 1
                cos_inc = np.ones(np.sum(mask_sol))

            elif tipus == "1-eix-horitzontal":
                # Cas 2: Beta fix, Gamma_p = Gamma_sol
                p_sol = p_unit[:, mask_sol]
                arg_sol = arg[mask_sol]
                e_est = np.array([-np.cos(arg_sol), -np.cos(delta_rad)*np.sin(arg_sol), np.sin(delta_rad)*np.sin(arg_sol)])
                e_nord = np.cross(p_sol.T, e_est.T).T

                # Calculem l'azimut del sol (gamma_s)
                d_E = np.dot(d_unit, e_est)
                d_N = np.dot(d_unit, e_nord)
                gamma_s = np.arctan2(d_E, -d_N)

                # Normal de la placa amb g_p = gamma_s
                b_rad = np.radians(b_fix)
                n_placa = (np.cos(b_rad) * p_sol +
                           np.sin(b_rad) * (-np.cos(gamma_s) * e_nord + np.sin(gamma_s) * e_est))
                cos_inc = np.einsum('i,ij->j', d_unit, n_placa)

            else: # Fixa 36° Sud
                p_sol = p_unit[:, mask_sol]
                arg_sol = arg[mask_sol]
                e_est = np.array([-np.cos(arg_sol), -np.cos(delta_rad)*np.sin(arg_sol), np.sin(delta_rad)*np.sin(arg_sol)])
                e_nord = np.cross(p_sol.T, e_est.T).T
                b_rad, g_rad = np.radians(b_fix), 0
                n_placa = (np.cos(b_rad) * p_sol +
                           np.sin(b_rad) * (-np.cos(g_rad) * e_nord + np.sin(g_rad) * e_est))
                cos_inc = np.einsum('i,ij->j', d_unit, n_placa)

            pot_inst = area * eficiencia * I_0 * np.maximum(0, cos_inc)
            potencies[mask_sol] = np.minimum(p_max, pot_inst)

        energia_dia_Wh = np.trapz(potencies, x=hores_dia)
        mitjanes_diaries.append(energia_dia_Wh / 24)

    return np.array(mitjanes_diaries)

# =============================================================================
# GENERACIÓ DE LA GRÀFICA
# =============================================================================
plt.figure(figsize=(15, 6))

# 1. Seguidor de 2 eixos (Cosinus = 1)
m_2e = simulacio_seguidors("2-eixos")
e_2e = np.trapz(m_2e, dx=1) * 24 / 1000
plt.plot(dies_any, m_2e, label=f"Seguidor 2 eixos (Cos θ = 1) [E={e_2e:.1f} kWh]", color='black')

# 2. Seguidor 1 eix (Beta=36, Gamma_p = Gamma_s)
m_1e = simulacio_seguidors("1-eix-horitzontal", b_fix=36)
e_1e = np.trapz(m_1e, dx=1) * 24 / 1000
plt.plot(dies_any, m_1e, label=f"Seguidor horitzontal (β=36°, γp=γs) [E={e_1e:.1f} kWh]", color='blue')

# 3. Referència: Fixa 36° Sud
m_fix = simulacio_seguidors("fix", b_fix=36)
e_fix = np.trapz(m_fix, dx=1) * 24 / 1000
plt.plot(dies_any, m_fix, label=f"Placa Fixa (β=36°, Sud) [E={e_fix:.1f} kWh]", color='green', alpha=0.6)

plt.tick_params(
    axis='both',
    direction='in',       # Guionets cap a dins
    top=True,             # Guionets al costat superior
    right=True,           # Guionets al costat dret
    which='both',
    labelsize=11
)

plt.xlim(0,400)
plt.ylim(0,300)


plt.xlabel("Dia de l'any", fontsize=12)
plt.ylabel("Potència Mitjana Diària (W)", fontsize=12)
plt.legend(fontsize=12)
plt.grid(alpha=0.3)
plt.show()