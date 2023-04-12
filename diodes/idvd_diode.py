# !/usr/bin/env python3
# -*- coding_ctrl: utf-8 -*-

# ------------------------------------- #
# Created by: JOaO PEDRO PETERS BARBOSA #
#           & MATHEUS SENE PAULO        #
#                                       #
# email: joao.peters@engenharia.ufjf.br #
#     & matheus.sene@engenharia.ufjf.br #
#                                       #
# Date: Jul/2021                        #
# ------------------------------------- #


"""
Disciplina [210081] - Tecnicas de Simulacao de Conversores Estaticos

Circuito de Cockcroft-Walton

Aplicacao da modelagem PWL no diodo, e metodo de integracao trapezoidal.

Baseado em material disponibilizado pelo professor.

Prof.: Pedro Gomes Barbosa

"""


import matplotlib.pyplot as plt
import numpy as np


Is = 1e-9
Vt = 25e-3

id, vd = 0, 0

n = np.arange(0, 2+10e-6, 10e-6)
outv = np.zeros(n.shape[0])
outi = np.zeros(n.shape[0])

for nn in range(n.shape[0]):
    id = Is * (np.exp(vd/Vt) - 1)
    vd += 10e-6

    outi[nn] = id
    outv[nn] = vd

print(vd)

plt.figure()
ax = plt.gca()
# plt.title("Curva Característica IxV Diodo")
plt.plot(outv+0.1, outi, label="Valor Real")
plt.xlim([0, 1.])
plt.ylim([-10, 100])

plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.plot([0, 0.7], [0, 1/1e6], label="Aproximação")
plt.plot([0.7, 0.8], [0, 1/5e-3], color='#F97306')
plt.legend(frameon=False)

plt.ylabel("Corrente $I_d$ [A]", fontsize=12)
plt.xlabel("Tensão $V_d$ [V]", fontsize=12)

plt.show()

