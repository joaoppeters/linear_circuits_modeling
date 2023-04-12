# !/usr/bin/env python3
# -*- coding_ctrl: utf-8 -*-

# ------------------------------------- #
# Created by: JOaO PEDRO PETERS BARBOSA #
#                                       #
# email: joao.peters@engenharia.ufjf.br #
# ------------------------------------- #

"""
Disciplina [210081] - Tecnicas de Simulacao de Conversores Estaticos

Resultados de circuito simples contendo diodo.

Ondas segundo equacao de Shockley.

Baseado em material disponibilizado pelo professor pgomes.
"""

import numpy as np
import matplotlib.pyplot as plt

Io = 1e-9

Vt = 25e-3

f = 1e3
T = 1 / f
w = 2 * np.pi * f

deltaT = T / 100
t = np.arange(0, 2 * T + deltaT, deltaT)

vs = 127 * np.sqrt(2) * np.sin(w * t)

R = 1e3
nptos = t.shape[0]

tol = 1e-6
vdini = 0.8
idk = Io * (np.exp(vdini / Vt) - 1)

vd = []
id = []

for i in range(nptos):
    errok = 1
    x = 0

    while (errok >= tol):
        D_vdini = -(R * Io / Vt) * np.exp(vdini / Vt) - 1
        F_vdini = vs[i] - R * idk - vdini
        vdk = vdini - F_vdini / D_vdini
        idk = Io * (np.exp(vdk / Vt) - 1)
        errok = np.abs(vdk - vdini)
        vdini = vdk

    vd.append(vdk)
    id.append(idk)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1)

ax1.plot(t, vs)
ax2.plot(t, vd)
ax3.plot(t, id)

plt.show()
