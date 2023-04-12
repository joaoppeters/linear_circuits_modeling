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
Jc3_ant2 = 0
Jc4_ant2 = 0
Jc5_ant2 = 0
Jc6_ant2 = 0

# -------------------- #
# Fonte de Alimentacao #
vS = 127  # (Vrms)
f = 60  # (Hertz)
w = 2 * np.pi * f  # (rad/s)
T = 1 / f  # (seg)

# ----------------------- #
# Parametros de Simulacao #
t0 = 0  # valor inicial
tn = 48 * T  # valor final
deltaT = 10e-6  # passo de simulacao

tempo = np.arange(t0, tn+deltaT, deltaT)  # tempo de simulacao
nptos = tempo.shape[0]  # total de pontos simulado

t_ant = 0  # tempo anterior (t - deltaT)

tol = 1e-3  # tolerancia de convergencia

# ----------------------------------------------- #
# Matriz Nodal Modificada (MNM) e Vetores Solucao #
MNM_var = np.zeros((8, 8))
MNM_cte = np.zeros((8, 8))

vi_k = np.zeros((8, 1))  # vetor solucao: v1, v2, v3, v4, v5, v6, v7 & iS
JE_k = np.zeros((8, 1))  # Contribuicoes fonte de corrente

vi_ant = vi_k  # vetor solucao anterior

# ------------------------------------ #
# Discretizacao dos elementos passivos #
Gc = (2 * C) / deltaT

# --------------------------------- #
# MNM constante: diagonal principal #
MNM_cte[0][0] = + Gc       # c1: 1->2
MNM_cte[1][1] = + Gc + Gc  # c1: 1->2 & c3: 4->2
MNM_cte[2][2] = + Gc + Gc  # c2: 3->g & c4: 5->3
MNM_cte[3][3] = + Gc + Gc  # c3: 4->2 & c5: 6->4
MNM_cte[4][4] = + Gc + Gc  # c4: 5->3 & c6: 7->5
MNM_cte[5][5] = + Gc       # c5: 6->4
MNM_cte[6][6] = + Gc       # c6: 7->5

# ------------------------------- #
# MNM constante: fora da diagonal #
MNM_cte[0][1] = - Gc  # c1: 1->2
MNM_cte[1][0] = - Gc  # :
MNM_cte[1][3] = - Gc  # c3: 4->2
MNM_cte[3][1] = - Gc  # :
MNM_cte[2][4] = - Gc  # c4: 5->3
MNM_cte[4][2] = - Gc  # :
MNM_cte[3][5] = - Gc  # c5: 6->4
MNM_cte[5][3] = - Gc  # :
MNM_cte[4][6] = - Gc  # c6: 7->5
MNM_cte[6][4] = - Gc  # :
MNM_cte[7][0] = + 1   # vS: 1->g
MNM_cte[0][7] = + 1   # :

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
out_v4 = np.zeros((1, nptos))
out_v5 = np.zeros((1, nptos))
out_v6 = np.zeros((1, nptos))
out_v7 = np.zeros((1, nptos))
out_iS = np.zeros((1, nptos))

# ----------------- #
# Loop de simulacao #
for t in range(1, nptos):

    #------------------------- #
    # atualizacao de variaveis #
    t_k = (deltaT * t)
    cont_d = 0
    errok = np.ones((8, 1))  # armazena erros

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

    # ---------------------------------------------------- #
    # Analise da tensao e fonte de corrente do capacitor 3 #
    vc3_ant = + vi_ant[3][0] - vi_ant[1][0]  # v4 - v2
    ic3_ant =  Gc * vc3_ant - Jc3_ant2
    Jc3_ant = + (Gc * vc3_ant) + ic3_ant  # ic + ic_ant -> tpz
    Jc3_ant2 = Jc3_ant

    # ---------------------------------------------------- #
    # Analise da tensao e fonte de corrente do capacitor 4 #
    vc4_ant = + vi_ant[4][0] - vi_ant[2][0]  # v5 - v3
    ic4_ant =  Gc * vc4_ant - Jc4_ant2
    Jc4_ant = + (Gc * vc4_ant) + ic4_ant  # ic + ic_ant -> tpz
    Jc4_ant2 = Jc4_ant

    # ---------------------------------------------------- #
    # Analise da tensao e fonte de corrente do capacitor 5 #
    vc5_ant = + vi_ant[5][0] - vi_ant[3][0]  # v6 - v4
    ic5_ant =  Gc * vc5_ant - Jc5_ant2
    Jc5_ant = + (Gc * vc5_ant) + ic5_ant  # ic + ic_ant -> tpz
    Jc5_ant2 = Jc5_ant

    # ---------------------------------------------------- #
    # Analise da tensao e fonte de corrente do capacitor 6 #
    vc6_ant = + vi_ant[6][0] - vi_ant[4][0]  # v7 - v5
    ic6_ant =  Gc * vc6_ant - Jc6_ant2
    Jc6_ant = + (Gc * vc6_ant) + ic6_ant  # ic + ic_ant -> tpz
    Jc6_ant2 = Jc6_ant

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
        vd2_ant = + vi_ant[1][0] - vi_ant[2][0]  # v2 - v3
        Gd2 = b + (cj * np.sign(vd2_ant - vD))
        id2 = a + (b * vd2_ant) + (cj * np.abs(vd2_ant - vD))
        Jd2 = id2 - (vd2_ant * Gd2)

        # ------------------------------------------------ #
        # Analise de tensao e fonte de corrente do diodo 3 #
        vd3_ant = + vi_ant[2][0] - vi_ant[3][0]  # v3 - v4
        Gd3 = b + (cj * np.sign(vd3_ant - vD))
        id3 = a + (b * vd3_ant) + (cj * np.abs(vd3_ant - vD))
        Jd3 = id3 - (vd3_ant * Gd3)

        # ------------------------------------------------ #
        # Analise de tensao e fonte de corrente do diodo 4 #
        vd4_ant = + vi_ant[3][0] - vi_ant[4][0]  # v4 - v5
        Gd4 = b + (cj * np.sign(vd4_ant - vD))
        id4 = a + (b * vd4_ant) + (cj * np.abs(vd4_ant - vD))
        Jd4 = id4 - (vd4_ant * Gd4)

        # ------------------------------------------------ #
        # Analise de tensao e fonte de corrente do diodo 5 #
        vd5_ant = + vi_ant[4][0] - vi_ant[5][0]  # v5 - v6
        Gd5 = b + (cj * np.sign(vd5_ant - vD))
        id5 = a + (b * vd5_ant) + (cj * np.abs(vd5_ant - vD))
        Jd5 = id5 - (vd5_ant * Gd5)

        # ------------------------------------------------ #
        # Analise de tensao e fonte de corrente do diodo 6 #
        vd6_ant = + vi_ant[5][0] - vi_ant[6][0]  # v6 - v7
        Gd6 = b + (cj * np.sign(vd6_ant - vD))
        id6 = a + (b * vd6_ant) + (cj * np.abs(vd6_ant - vD))
        Jd6 = id6 - (vd6_ant * Gd6)

        # -------------------------------- #
        # MNM variavel: diagonal principal #
        MNM_var[1][1] = + Gd1 + Gd2  # d1: 2->g & d2: 2->3
        MNM_var[2][2] = + Gd2 + Gd3  # d2: 2->3 & d3: 3->4
        MNM_var[3][3] = + Gd3 + Gd4  # d3: 3->4 & d4: 4->5
        MNM_var[4][4] = + Gd4 + Gd5  # d4: 4->5 & d5: 5->6
        MNM_var[5][5] = + Gd5 + Gd6  # d5: 5->6 & d6: 6->7
        MNM_var[6][6] = + Gd6        # d6: 6->7

        # ------------------------------ #
        # MNM variavel: fora da diagonal #
        MNM_var[1][2] = - Gd2  # d2: 2->3
        MNM_var[2][1] = - Gd2  # :
        MNM_var[2][3] = - Gd3  # d3: 3->4
        MNM_var[3][2] = - Gd3  # :
        MNM_var[3][4] = - Gd4  # d4: 4->5
        MNM_var[4][3] = - Gd4  # :
        MNM_var[4][5] = - Gd5  # d5: 5->6
        MNM_var[5][4] = - Gd5  # :
        MNM_var[5][6] = - Gd6  # d6: 6->5
        MNM_var[6][5] = - Gd6  # :

        invMNM = np.linalg.inv(MNM_cte + MNM_var)  # MNM invertida

        # ----------------------------- #
        # Vetor de Fontes Independentes #
        JE_k[0][0] = + Jc1_ant
        JE_k[1][0] = - Jc1_ant - Jc3_ant + Jd1 - Jd2
        JE_k[2][0] = + Jc2_ant - Jc4_ant + Jd2 - Jd3
        JE_k[3][0] = + Jc3_ant - Jc5_ant + Jd3 - Jd4
        JE_k[4][0] = + Jc4_ant - Jc6_ant + Jd4 - Jd5
        JE_k[5][0] = + Jc5_ant + Jd5 - Jd6
        JE_k[6][0] = + Jc6_ant + Jd6
        JE_k[7][0] = + vs_k

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
    out_v4[0][t] = vi_k[3][0]
    out_v5[0][t] = vi_k[4][0]
    out_v6[0][t] = vi_k[5][0]
    out_v7[0][t] = vi_k[6][0]
    out_iS[0][t] = vi_k[7][0]

# --------------------------------- #
# valor maximo de tensao & corrente #
print(np.max(out_v7[0][:]), np.max(out_iS[0][:]))

# -------- #
# FIGURA 1 #
plt.figure()
plt.plot(tempo, out_v1[0][:], label="Tensao Fonte")
plt.plot(tempo, out_v2[0][:], label="Tensao N2")
plt.plot(tempo, out_v3[0][:], label="Tensao N3")
plt.plot(tempo, out_v4[0][:], label="Tensao N4")
plt.plot(tempo, out_v5[0][:], label="Tensao N5")
plt.plot(tempo, out_v6[0][:], label="Tensao N6")
plt.plot(tempo, out_v7[0][:], label="Tensao N7")
plt.plot(tempo, out_iS[0][:], label="Corrente Fonte")

plt.xlabel("Tempo de Simulacao (seg)")
plt.ylabel("Amplitude")

plt.legend(frameon=False)

# -------- #
# FIGURA 2 #
plt.figure()
plt.plot(tempo, out_iS[0][:], label="Corrente Fonte")

plt.legend(frameon=False)
plt.show()