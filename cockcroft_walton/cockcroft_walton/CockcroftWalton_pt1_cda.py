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
import os
import sys

# ------------------------------ #
# Elementos Passivos do Circuito #
C = 500e-6  # (Farad)

# -------------------------------------- #
# Fonte de Corrente dos Capacitores: TPZ #
Jc1_ant2 = 0
Jc2_ant2 = 0

# -------------------- #
# Fonte de Alimentacao #
vS = 127  # (Vrms)
f = 60  # (Hertz)
w = 2 * np.pi * f  # (rad/s)
T = 1 / f  # (seg)

# ----------------------- #
# Parametros de Simulacao #
t0 = 0  # valor inicial
tn = 12 * T  # valor final
deltaT = 100e-6  # passo de simulacao

tempo = np.arange(t0, tn+deltaT, deltaT)  # tempo de simulacao
nptos = tempo.shape[0]  # total de pontos simulado

t_ant = 0  # tempo anterior (t - deltaT)

tol = 1e-3  # tolerancia de convergencia

# ----------------------------------------------- #
# Matriz Nodal Modificada (MNM) e Vetores Solucao #
# MNM = np.zeros((3, 3))
MNM_var = np.zeros((4, 4))
MNM_cte = np.zeros((4, 4))

vi_k = np.zeros((4, 1))  # vetor solucao: v1, v2, v3 & iS
JE_k = np.zeros((4, 1))  # Contribuicoes fonte de corrente

vi_ant = vi_k  # vetor solucao anterior

# ------------------------------------ #
# Discretizacao dos elementos passivos #
Gc = (2 * C) / deltaT

# --------------------------------- #
# MNM constante: diagonal principal #
MNM_cte[0][0] = + Gc  # c1: 1->2
MNM_cte[1][1] = + Gc  # :
MNM_cte[2][2] = + Gc  # c2: 3->g

# ------------------------------- #
# MNM constante: fora da diagonal #
MNM_cte[1][0] = - Gc  # c1: 1->2
MNM_cte[0][1] = - Gc  # :
MNM_cte[3][0] = + 1   # vS: 1->g
MNM_cte[0][3] = + 1   # :

# ---------------------- #
# Modelagem PWL do diodo #
Rd0 = 1e6
Rd1 = 10e-3

Gd0 = 1 / Rd0
Gd1 = 1 / Rd1

vD = 0.7
i0 = 0

b = (Gd1 + Gd0) / 2
cj = (Gd1 - Gd0) / 2
a = i0 - (cj * np.abs(vD))

# ------------------------------- #
# Variaveis de saida da simulacao #
out_v1 = np.zeros((1, nptos))
out_v2 = np.zeros((1, nptos))
out_v3 = np.zeros((1, nptos))
out_iS = np.zeros((1, nptos))

# ----------------- #
# Loop de simulacao #
for t in range(1, nptos):

    #------------------------- #
    # atualizacao de variaveis #
    t_k = (deltaT * t)
    cont_d = 0
    errok = np.ones((4, 1))  # armazena erros

    # --------------- #
    # fonte de tensao #
    vs_k = np.sqrt(2) * vS * np.sin((w * t_k) + np.pi)  # defasagem 180g

    # ---------------------------------------------------- #
    # Analise da tensao e fonte de corrente do capacitor 1 #
    vc1_ant = + vi_ant[0][0] - vi_ant[1][0]  # v1 - v2
    ic1_ant =  Gc * vc1_ant - Jc1_ant2
    Jc1_ant = + (Gc * vc1_ant) + ic1_ant  # ic + ic_ant -> tpz
    Jc1_ant2 = Jc1_ant

    # ---------------------------------------------------- #
    # Analise da tensao e fonte de corrente do capacitor 2 #
    vc2_ant = + vi_ant[2][0]  # v3 - g
    ic2_ant =  Gc * vc2_ant - Jc2_ant2
    Jc2_ant = + (Gc * vc2_ant) + ic2_ant  # ic + ic_ant -> tpz
    Jc2_ant2 = Jc2_ant

    # ------------------- #
    # Modelagem PWL diodo #
    while np.max(np.abs(errok)) > tol:
        # print(t, cont_d)

        # ------------------------------------------------ #
        # Analise de tensao e fonte de corrente do diodo 1 #
        vd1_ant = - vi_ant[1][0]  # g - v2
        Gd1 = b + (cj * np.sign(vd1_ant - vD))
        id1 = a + (b * vd1_ant) + (cj * np.abs(vd1_ant - vD))
        Jd1 = id1 - (vd1_ant * Gd1)

        # ------------------------------------------------ #
        # Analise de tensao e fonte de corrente do diodo 2 #
        vd2_ant = vi_ant[1][0] - vi_ant[2][0]  # v2 - v3
        Gd2 = b + (cj * np.sign(vd2_ant - vD))
        id2 = a + (b * vd2_ant) + (cj * np.abs(vd2_ant - vD))
        Jd2 = id2 - (vd2_ant * Gd2)

        # -------------------------------- #
        # MNM variavel: diagonal principal #
        MNM_var[1][1] = + Gd1 + Gd2  # d1: 2->g & d2: 2->3
        MNM_var[2][2] = + Gd2        # d2: 2->3

        # ------------------------------ #
        # MNM variavel: fora da diagonal #
        MNM_var[1][2] = - Gd2        # d2: 2->3
        MNM_var[2][1] = - Gd2        # :

        invMNM = np.linalg.inv(MNM_cte + MNM_var)  # MNM invertida

        # ----------------------------- #
        # Vetor de Fontes Independentes #
        JE_k[0][0] = + Jc1_ant
        JE_k[1][0] = + Jd1 - Jd2 - Jc1_ant
        JE_k[2][0] = + Jd2 + Jc2_ant
        JE_k[3][0] = + vs_k

        # -------------------------- #
        # Calculo das tensoes nodais #
        vi_k = np.dot(invMNM, JE_k)
        
        errok = (vi_k - vi_ant)
        vi_ant = vi_k
        
        cont_d += 1
        if cont_d > 100:
            print("erro, iteracao {}".format(t+1))
            sys.exit()

    # ------------------------ #
    # Armazenamento de valores #
    out_v1[0][t] = vi_k[0][0]
    out_v2[0][t] = vi_k[1][0]
    out_v3[0][t] = vi_k[2][0]
    out_iS[0][t] = vi_k[3][0]
        
# --------------------------------- #
# valor maximo de tensao & corrente #
print(np.max(out_v3[0][:]), np.max(out_iS[0][:]))

# -------- #
# FIGURA 1 #
plt.figure()
plt.plot(tempo, out_v1[0][:], label="Tensao Fonte")
plt.plot(tempo, out_v2[0][:], label="Tensao N2")
plt.plot(tempo, out_v3[0][:], label="Tensao N3")
plt.plot(tempo, out_iS[0][:], label="Corrente Fonte")

plt.xlabel("Tempo de Simulacao (seg)")
plt.ylabel("Amplitude")
plt.legend(frameon=True, facecolor="white", edgecolor="white")

# -------- #
# FIGURA 2 #
plt.figure()
plt.plot(tempo, out_iS[0][:], label="Corrente Fonte", color="grey")

plt.xlabel("Tempo de Simulacao (seg)")
plt.ylabel("Amplitude [A]")
plt.legend(frameon=True, facecolor="white", edgecolor="white")
plt.grid()
plt.show()