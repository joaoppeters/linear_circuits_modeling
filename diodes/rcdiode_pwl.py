# !/usr/bin/env python3
# -*- coding_ctrl: utf-8 -*-

# ------------------------------------- #
# Created by: JOaO PEDRO PETERS BARBOSA #
#                                       #
# email: joao.peters@engenharia.ufjf.br #
# ------------------------------------- #


"""
Disciplina [210081] - Tecnicas de Simulacao de Conversores Estaticos

Circuito RLC: Indutor substituido por diodo.

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
R = 12  # (Ohms)
C = 500e-6  # (Farad)

# -------------------- #
# Fonte de Alimentacao #
vS = 127  # (Vrms)
f = 60  # (Hertz)
w = 2 * np.pi * f  # (rad/s)
T = 1 / f  # (seg)

# ----------------------- #
# Parametros de Simulacao #
t0 = 0  # valor inicial
tn = 3 * T  # valor final
deltaT = 10e-6  # passo de simulacao

tempo = np.arange(t0, tn+deltaT, deltaT)  # tempo de simulacao
nptos = tempo.shape[0]  # total de pontos simulado

t_ant = 0  # tempo anterior (t - deltaT)

tol = 1e-3  # tolerancia de convergencia

# ----------------------------------------------- #
# Matriz Nodal Modificada (MNM) e Vetores Solucao #
# MNM = np.zeros((3, 3))
MNM_var = np.zeros((3, 3))
MNM_cte = np.zeros((3, 3))

vi_k = np.zeros((3, 1))  # vetor solucao
JE_k = np.zeros((3, 1))  # Contribuicoes fonte de corrente

vi_ant = vi_k  # vetor solucao anterior

# ------------------------------------ #
# Discretizacao dos elementos passivos #
Gc = (2 * C) / deltaT
Gr = 1 / R

# ------------- #
# MNM constante #
MNM_cte[1][1] = Gc + Gr
MNM_cte[2][0] = 1
MNM_cte[0][2] = 1

#MNM_cte = np.array([[ 0.,     0.  ,+1.],
#                    [ 0.,  (Gc+Gr), 0.],
#                    [+1.,     0.  , 0.]])


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
out_iS = np.zeros((1, nptos))

# ----------------- #
# Loop de simulacao #
for t in range(1, nptos):

    #------------------------- #
    # atualizacao de variaveis #
    t_k = (deltaT * t)
    cont_d = 0
    errok = np.ones((3, 1))  # armazena erros

    # --------------- #
    # fonte de tensao #
    vs_k = np.sqrt(2) * vS * np.sin(w * t_k)

    # -------------------------------------------------- #
    # Analise da tensao e fonte de corrente do capacitor #
    vc_ant = + vi_ant[1][0]
    ic_ant = - vi_ant[2][0] - (Gr * vc_ant)
    Jc_ant = + (Gc * vc_ant) + ic_ant
    
    # --------- #
    # PWL diodo #
    while np.max(np.abs(errok)) > tol:
        # print(t, cont_d)

        # ---------------------------------------------- #
        # Analise de tensao e fonte de corrente do diodo #
        vd_ant = vi_ant[0][0] - vi_ant[1][0]  # v1 - v2
        Gd = b + (cj * np.sign(vd_ant - vD))
        id = a + (b * vd_ant) + (cj * np.abs(vd_ant - vD))
        Jd = id - (vd_ant * Gd)

        # -------------------------------- #
        # Matriz Nodal Modificada variavel #
        MNM_var[0][0] = +Gd
        MNM_var[0][1] = -Gd
        MNM_var[1][0] = -Gd
        MNM_var[1][1] = +Gd
        
        invMNM = np.linalg.inv(MNM_cte + MNM_var)  # MNM invertida

        # ----------------------------- #
        # Vetor de Fontes Independentes #
        JE_k[0][0] = - Jd
        JE_k[1][0] = + Jd + Jc_ant
        JE_k[2][0] = + vs_k

        # -------------------------- #
        # Calculo das tensoes nodais #
        vi_k = np.dot(invMNM, JE_k)
        
        errok = (vi_k - vi_ant)
        vi_ant = vi_k
        
        cont_d += 1
        if cont_d > 100:
            print("erro, iteracao {}".format(t+1))
            sys.exit()

    # t_ant = t_k

    out_v1[0][t] = vi_k[0][0]
    out_v2[0][t] = vi_k[1][0]
    out_iS[0][t] = vi_k[2][0]
        
    
plt.figure()
plt.plot(tempo, out_v1[0][:], label="Tensao Fonte")
plt.plot(tempo, out_v2[0][:], label="Tensao Capacitor")
plt.plot(tempo, out_iS[0][:], label="Corrente Fonte")

plt.xlabel("Tempo de Simulacao (seg)")
plt.ylabel("Amplitude")

plt.legend(frameon=False)
plt.show()
