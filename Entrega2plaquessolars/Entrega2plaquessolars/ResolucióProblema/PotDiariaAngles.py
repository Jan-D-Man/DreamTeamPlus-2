import numpy as np
import matplotlib.pyplot as plt
from Definicions import *



day_of_year = 100 # Dia de primavera random
LSTs = np.linspace(0, 24, 500)

# Paràmetres placa solar
beta_fixed = np.pi/4  # 45 graus
angles_or_graus = [-90, -30, -15, 0, 15, 30, 90]
cmap = plt.get_cmap('tab10')

plt.figure(figsize=(10, 5))

for i, g_graus in enumerate(angles_or_graus):
    g_rad = np.radians(g_graus)
    pot_dia = []

    # Càlculs del dia
    theta_dia = (2 * np.pi * (day_of_year - 172) / 365) - (np.pi / 2)
    t_central = (theta_dia / omega_T) + (np.pi / (2 * omega_T))

    for lst in LSTs:
        t = t_central + (lst - 12)
        arg = omega_T * t

        # Vectors Rodrigues
        p_vec = np.array([
            -RT * np.sin(mu) * np.sin(arg),
            RT * (np.cos(mu)*np.sin(delta_rad) + np.sin(mu)*np.cos(delta_rad)*np.cos(arg)),
            RT * (np.cos(mu)*np.cos(delta_rad) - np.sin(mu)*np.sin(delta_rad)*np.cos(arg))
        ])
        p_unit = p_vec / RT

        # Sol
        s_vec = np.array([RS * np.cos(theta_dia), RS * np.sin(theta_dia), 0])
        d_vec = -p_vec - s_vec
        d_unit = d_vec / np.linalg.norm(d_vec)

        # Elevació alpha
        alpha = np.arcsin(np.clip(np.dot(d_unit, p_unit), -1, 1))

        pot = 0
        if alpha > 0:
            e_est = np.array([-np.cos(arg), -np.cos(delta_rad)*np.sin(arg), np.sin(delta_rad)*np.sin(arg)])
            e_nord = np.cross(p_unit, e_est)

            # La teva n_placa generalitzada
            n_placa = (np.cos(beta_fixed) * p_unit +
                       np.sin(beta_fixed) * (-np.cos(g_rad) * e_nord + np.sin(g_rad) * e_est))

            cos_inc = np.dot(d_unit, n_placa)
            if cos_inc > 0:
                pot = min(p_max, area * eficiencia * I_0 * cos_inc)

        pot_dia.append(pot)

    energia = np.trapz(pot_dia, x=LSTs)
    plt.plot(LSTs, pot_dia, color=cmap(i), label=f"{g_graus}° (E={energia:.0f} Wh)")

plt.tick_params(
    axis='both',
    direction='in',       # Guionets cap a dins
    top=True,             # Guionets al costat superior
    right=True,           # Guionets al costat dret
    which='both',
    labelsize=11
)

plt.xlim(0,25)
plt.ylim(0,450)


plt.xlabel("Hora solar (LST)"); plt.ylabel("Potència (W)")
plt.legend(title="Azimut placa", bbox_to_anchor=(1.05, 1), fontsize=12, title_fontsize=13)
plt.grid(True, alpha=0.3); plt.tight_layout()

# comparem amb diferents angles beta fixant l'azimut a 0° (Sud)
betas_graus = [0, 15, 30, 45, 60, 90]
g_fixed = 0 # Sud pur

plt.figure(figsize=(10, 5))

for i, b_graus in enumerate(betas_graus):
    b_rad = np.radians(b_graus)
    pot_dia = []

    # Càlculs del dia (es repeteixen per claredat)
    theta_dia = (2 * np.pi * (day_of_year - 172) / 365) - (np.pi / 2)
    t_central = (theta_dia / omega_T) + (np.pi / (2 * omega_T))

    # Càlcul de la potència al llarg del dia
    for lst in LSTs:
        t = t_central + (lst - 12)
        arg = omega_T * t
        # Vectors Rodrigues
        p_vec = np.array([
            -RT * np.sin(mu) * np.sin(arg),
            RT * (np.cos(mu)*np.sin(delta_rad) + np.sin(mu)*np.cos(delta_rad)*np.cos(arg)),
            RT * (np.cos(mu)*np.cos(delta_rad) - np.sin(mu)*np.sin(delta_rad)*np.cos(arg))
        ])
        p_unit = p_vec / RT
        s_vec = np.array([RS * np.cos(theta_dia), RS * np.sin(theta_dia), 0])
        d_unit = (-p_vec - s_vec) / np.linalg.norm(-p_vec - s_vec)

        alpha = np.arcsin(np.clip(np.dot(d_unit, p_unit), -1, 1))
        # Càlcul de la potència si és de dia
        pot = 0
        if alpha > 0:
            e_est = np.array([-np.cos(arg), -np.cos(delta_rad)*np.sin(arg), np.sin(delta_rad)*np.sin(arg)])
            e_nord = np.cross(p_unit, e_est)
            n_placa = (np.cos(b_rad) * p_unit +
                       np.sin(b_rad) * (-np.cos(g_fixed) * e_nord + np.sin(g_fixed) * e_est))

            cos_inc = np.dot(d_unit, n_placa)
            if cos_inc > 0:
                pot = min(p_max, area * eficiencia * I_0 * cos_inc)
        pot_dia.append(pot)

    energia = np.trapz(pot_dia, x=LSTs) # Energia del dia
    plt.plot(LSTs, pot_dia, color=cmap(i), label=f"β={b_graus}° (E={energia:.0f} Wh)") 
plt.tick_params(
    axis='both',
    direction='in',       # Guionets cap a dins
    top=True,             # Guionets al costat superior
    right=True,           # Guionets al costat dret
    which='both',
    labelsize=11
)

plt.xlim(0,25)
plt.ylim(0,450)
plt.xlabel("Hora solar (LST)"); plt.ylabel("Potència (W)")
plt.legend(title="Inclinació β", bbox_to_anchor=(1.05, 1), fontsize=12, title_fontsize=13)
plt.grid(True, alpha=0.3); plt.tight_layout()
plt.show()