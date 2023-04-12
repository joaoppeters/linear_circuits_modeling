#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
Simulacao digital de um circuito multiplicador (Cockroft-Walton Multiplier) excitado por uma fonte senoidal de 60 Hz

Created on Sun Jul  15 21:22:53 2021

@author: pgomes
"""

import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------ #
# Definicao dos parametros do circuito #
# ------------------------------------ #
C = 500.0e-6
Jc1_ant = 0.0
Jc2_ant = 0.0

# --------------------------------------- #
# Fonte senoidal de excitacao do circuito #
# --------------------------------------- #
f = 60.0
T = 1 / f
Vs = 127.
omega = 2 * np.pi * f

# ------------------------------------ #
# Definicao de parametros da simulacao #
# ------------------------------------ #
t_ini = 0.  # tempo de inicio da simulacao
t_end = 0.2  # tempo de final de simulacao
deltaT = 100.0e-6  # passo de simulacao/integracao = deltaT

tempo = np.arange(t_ini, t_end + deltaT, deltaT)  # vetor de tempo da simulacao
nptos = tempo.shape[-1]  # numero de pontos da simulacao

t_k = t_ini  # tempo atual
t_ant = t_ini  # tempo anterior (t-deltaT)

n = 0  # contador do tempo discreto
j = 0  # contador do número de iterações
tol = 1.0e-2  # tolerância

# ------------------------------------ #
# Parametros para implementação do CDA #
# ------------------------------------ #

h = deltaT
cda = 0
be = 0

# ------------------------------------ #
# Definicao da matriz nodal modificada #
# ------------------------------------ #
MNM_cte = np.zeros((4, 4))
MNM_var = np.zeros((4, 4))
vn = np.zeros((4, 1))
JE_k = np.zeros((4, 1))

vn_ant = vn

# ------------------------------------ #
# Discretizacao dos elementos passivos #
Gc = (2 * C) / deltaT

# ------------- #
# MNM constante #
MNM_cte[0][0] = + Gc  # c1: 1->2
MNM_cte[1][0] = - Gc  # :
MNM_cte[0][1] = - Gc  # :
MNM_cte[1][1] = + Gc  # :
MNM_cte[2][2] = + Gc  # c2: 3->gr
MNM_cte[3][0] = + 1   # vS: 1->gr
MNM_cte[0][3] = + 1   # :

# -----------------------------------------------------
# aproximacao pwl para o diodo
# -----------------------------------------------------

# # utilização da equação de Schokly para determinar as inclinações da curva do diodo
VT = 25.0e-3  # tensão térmica do diodo real
I0 = 1.0e-9  # corrente de fuga do diodo real

G0 = 1.0e-6
G1 = 2.0e+03
V1 = 0.0

b = (G0 + G1) / 2
Cj = (G1 - G0) / 2
a = I0 - Cj * np.abs(V1)

# ----------------------------------------------#
# Vetor com as condutâncias iniciais dos diodos #
# ----------------------------------------------#

nDiodos = 2
GD_k = np.zeros((nDiodos, 1))
GD_ant = np.zeros((nDiodos, 1))

# -------------------------------------------------------------------- #
# saida de variaveis ( tensoes nodais + corrente pela fonte de tensao) #
# -------------------------------------------------------------------- #
out1 = np.zeros((nptos))
out2 = np.zeros((nptos))
out3 = np.zeros((nptos))
out4 = np.zeros((nptos))

# ----------------- #
# Loop de simulacao #
# ----------------- #
while t_k < t_end:

    erro_k = np.ones((4, 1))  # vetor de erro da iteração k
    # -----------------------------#
    # Fonte de tensao independente #
    # -----------------------------#
    vs_k = np.sqrt(2) * Vs * np.sin(omega * t_k + np.pi)

    # -----------------------------------------------#
    # Cálculo das fontes de corrente dos capacitores #
    # -----------------------------------------------#

    # Capacitor C1
    vc1_ant = (vn_ant[0][0] - vn_ant[1][0])  # tensao do no 2 do passado
    ic1_ant = (Gc * vc1_ant - Jc1_ant) * (1 - be)
    Jc1_ant = Gc * vc1_ant + ic1_ant * (1 - be)

    # Capacitor C2
    vc2_ant = (vn_ant[2][0])  # tensao do no 2 do passado
    ic2_ant = (Gc * vc2_ant - Jc2_ant) * (1 - be)
    Jc2_ant = Gc * vc2_ant + ic2_ant * (1 - be)

    vnk_ = vn_ant

    while np.max(np.abs(erro_k)) >= tol:
        # -------------------------------------------------------#
        # Cálculo da condutância e fonte de corrente dos diodos  #
        # -------------------------------------------------------#

        # Contribuição do diodo D1
        vd1_ant = 0 - vn_ant[1][0]
        GD1k_ant = b + Cj * np.sign(vd1_ant - V1)
        ID1k_ant = a + b * vd1_ant + Cj * np.abs(vd1_ant - V1)
        Ieq1k_ant = ID1k_ant - vd1_ant * GD1k_ant

        # Contribuição do diodo D2
        vd2_ant = vn_ant[1][0] - vn_ant[2][0]
        GD2k_ant = b + Cj * np.sign(vd2_ant - V1)
        ID2k_ant = a + b * vd2_ant + Cj * np.abs(vd2_ant - V1)
        Ieq2k_ant = ID2k_ant - vd2_ant * GD2k_ant

        # -------------------------------- #
        # Matriz Nodal Modificada variavel #
        MNM_var[1][1] = + GD1k_ant + GD2k_ant  # d1: 2->gr || d2:2->3
        MNM_var[1][2] = - GD2k_ant  # :
        MNM_var[2][1] = - GD2k_ant  # :
        MNM_var[2][2] = + GD2k_ant  # :

        invMNM = np.linalg.inv(MNM_cte + MNM_var)  # MNM invertida

        # ----------------------------- #
        # Vetor de Fontes Independentes #
        JE_k[0][0] = + Ieq1k_ant
        JE_k[1][0] = + ID1k_ant - ID2k_ant - Ieq1k_ant
        JE_k[2][0] = + ID2k_ant + Ieq2k_ant
        JE_k[3][0] = + vs_k

        # -------------------------- #
        # Calculo das tensoes nodais #
        # -------------------------- #
        vnk = np.dot(invMNM, JE_k)

        erro_k = (vnk - vnk_)
        vnk_ = vnk

        j = j + 1
        print("iter = ", j)

    # ----------------------------------------- #
    # Atualizacao das tensoes nodais e contador #
    # ----------------------------------------- #
    GD_k[0][0] = GD1k_ant
    GD_k[1][0] = GD2k_ant

    if ((GD_k != GD_ant).any()):
        if be == 0:
            be = 1
            cda = 3
            h = deltaT / 2
            t_k = t_ant - h
            vnk_ = vn_ant

    if cda != 0:
        j = 0

        t_ant = t_k
        t_k = t_ant + h
        cda = cda - 1
        GD_ant[0][0] = GD_k[0][0]
        GD_ant[1][0] = GD_k[1][0]

    if cda == 0:
        be = 0

        print("t_n = ", n)

        h = deltaT

        n = n + 1
        t_ant = t_k
        t_k = t_ant + h  # avança o tempo da simulacao

        j = 0

        vn = vnk
        vn_ant = vnk

        GD_ant[0][0] = GD_k[0][0]
        GD_ant[1][0] = GD_k[1][0]
        # ------------------ #    
        # Saida de variaveis #
        # ------------------ #    
        out1[n] = vn[0][0]
        out2[n] = vn[1][0]
        out3[n] = vn[2][0]
        out4[n] = vn[3][0]

plt.figure()
plt.plot(tempo, out1, tempo, out2, tempo, out3)
plt.xlabel('tempo [s]')
plt.legend(['V_1', 'V_2'])

plt.figure()
plt.plot(tempo, out4)
plt.xlabel('tempo [s]')
plt.legend(['i_s'])
plt.show()

