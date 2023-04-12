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
# R = 12  # (Ohms)
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
tn = 0.2  # valor final
deltaT = 100e-6  # passo de simulacao

tempo1 = np.arange(t0, tn+deltaT, deltaT)  # tempo de simulacao
nptos = tempo1.shape[0]  # total de pontos simulado

t_ant = 0  # tempo anterior (t - deltaT)
tol = 1e-3  # tolerancia de convergencia

# ----------------------------------------------- #
# Matriz Nodal Modificada (MNM) e Vetores Solucao #
# MNM = np.zeros((3, 3))
MNM_var = np.zeros((4, 4))
MNM_cte = np.zeros((4, 4))

vi = np.zeros((4, 1))  # vetor solucao: v1, v2, v3 & iS
JE_k = np.zeros((4, 1))  # Contribuicoes fonte de corrente

vi_ant = vi  # vetor solucao anterior
vk_ant = vi

# ---------------------- #
# Modelagem PWL do diodo #
Rd0 = 1e6
Rd1 = 10e-3

Gd0 = 1 / Rd0
Gd1 = 1 / Rd1
Gc = (2 * C) / deltaT

vD = 0.7
i0 = 0

b = (Gd1 + Gd0) / 2
cj = (Gd1 - Gd0) / 2
a = i0 - (cj * np.abs(vD))

# ------------------------------- #
# Variaveis de saida da simulacao #
out_v1 = [0]
out_v2 = [0]
out_v3 = [0]
out_iS = [0]
Jc2_ant2 = 0 
Jc1_ant2 = 0 
id1_vet = [0]
id2_vet = [0]
flag = 'TPZ'
passos = 0 
tempo= [0]
t = 1
Gd1_vet=[0]
Gd2_vet=[0]

    # MNM constante 
MNM_cte[0][0] = + Gc  # c1: 1->2
MNM_cte[1][0] = - Gc  # :
MNM_cte[0][1] = - Gc  # :
MNM_cte[1][1] = + Gc  # :
MNM_cte[2][2] = + Gc  # c2: 3->gr
MNM_cte[3][0] = + 1   # vS: 1->gr
MNM_cte[0][3] = + 1   # :       
# ----------------- #
# Loop de simulacao #
while t < nptos:
    cont_d = 0
    errok = np.ones((4, 1))  # armazena erros

    vc1_ant = + vi_ant[0][0] - vi_ant[1][0]  # v1 - v2
    vc2_ant = + vi_ant[2][0]  

    if flag == 'TPZ':

        t_k = (deltaT * t)
        tempo.append(t_k)
        # --------------- #
        # fonte de tensao #
        vs_k = np.sqrt(2) * vS * np.sin((w * t_k) + np.pi)
     
        ic1_ant =  Gc * vc1_ant - Jc1_ant2
        Jc1_ant = + (Gc * vc1_ant) + ic1_ant  # ic + ic_ant -> tpz
        Jc1_ant2 = Jc1_ant

        ic2_ant =  Gc * vc2_ant - Jc2_ant2
        Jc2_ant = + (Gc * vc2_ant) + ic2_ant  # ic + ic_ant -> tpz
        Jc2_ant2 = Jc2_ant

    if flag == 'BE':

        t_k = t_k + (deltaT/2)
        tempo.append(t_k)
        # --------------- #
        # fonte de tensao #
        vs_k = np.sqrt(2) * vS * np.sin((w * t_k) + np.pi)
              
        Jc1_ant = + (Gc * vc1_ant) 
        Jc1_ant2 = Jc1_ant

        Jc2_ant = + (Gc * vc2_ant)   
        Jc2_ant2 = Jc2_ant

        passos -= 1          
        
    # ------------------- #
    # Modelagem PWL diodo #

    while np.max(np.abs(errok)) > tol:
        print(t, cont_d)

        # ------------------------------------------------ #
        # Analise de tensao e fonte de corrente do diodo 1 #
        vd1_ant = - vk_ant[1][0]  # -v2
        Gd1 = b + (cj * np.sign(vd1_ant - vD))
        
        id1 = a + (b * vd1_ant) + (cj * np.abs(vd1_ant - vD))       
        Jd1 = id1 - (vd1_ant * Gd1)
        # ------------------------------------------------ #
        # Analise de tensao e fonte de corrente do diodo 2 #
        vd2_ant = vk_ant[1][0] - vk_ant[2][0]  # v2 - v3
        Gd2 = b + (cj * np.sign(vd2_ant - vD))
        
        id2 = a + (b * vd2_ant) + (cj * np.abs(vd2_ant - vD))
        #id2_vet.append(id2)
        Jd2 = id2 - (vd2_ant * Gd2)     

        # -------------------------------- #
        # Matriz Nodal Modificada variavel #
        MNM_var[1][1] = + Gd1 + Gd2  # d1: 2->gr || d2:2->3
        MNM_var[1][2] = - Gd2        # :
        MNM_var[2][1] = - Gd2        # :
        MNM_var[2][2] = + Gd2        # :

        invMNM = np.linalg.inv(MNM_cte + MNM_var)  # MNM invertida

        # ----------------------------- #
        # Vetor de Fontes Independentes #
        JE_k[0][0] = + Jc1_ant
        JE_k[1][0] = + Jd1 - Jd2 - Jc1_ant
        JE_k[2][0] = + Jd2 + Jc2_ant
        JE_k[3][0] = + vs_k

        # -------------------------- #
        # Calculo das tensoes nodais #
        vi = np.dot(invMNM, JE_k)
        
        cont_d = cont_d + 1
        errok = (vi - vk_ant) #vk_ant roda no loop iterativo
        vk_ant = vi       

    Gd1_vet.append(Gd1)
    Gd2_vet.append(Gd2)
    #id1_vet.append(id1) 
    #id2_vet.append(id2)
    
    d = (vi[3][0] - vi_ant[3][0])/deltaT

    #if ( (Gd1_vet[t] != Gd1_vet[t-1] or Gd2_vet[t] != Gd2_vet[t-1]) and flag == 'TPZ'): 
    if ( np.abs(d) > 20e3 and flag == 'TPZ'):
        flag = 'BE'
        passos = 2
        nptos += 2
        tempo.pop(t)
        Gd1_vet.pop(t)
        Gd2_vet.pop(t)
        t = t - 1 
        vk_ant = vi_ant  
        t_k = t_k - deltaT
             
    else: 
        vi_ant = vi
        out_v1.append(vi[0][0])
        out_v2.append(vi[1][0])
        out_v3.append(vi[2][0])
        out_iS.append(vi[3][0])
        t += 1

    if passos == 0  and flag == 'BE': 
        flag='TPZ'
        #t = t - 1  


plt.figure()
plt.plot(tempo, out_v1, label="Tensao Fonte", color="red")
plt.plot(tempo, out_v2, label="Tensao N2", color="lightgreen")
plt.plot(tempo, out_v3, label="Tensao N3", color="darkgreen")
plt.plot(tempo, out_iS, label="Corrente Fonte", color="grey")

plt.xlabel("Tempo de Simulacao (seg)", fontsize=12)
plt.ylabel("Amplitude [V] & [A]", fontsize=12)
plt.legend(frameon=True, facecolor="white", edgecolor="white")
plt.grid()


# FIGURA 2 #
plt.figure()
plt.plot(tempo, out_iS, label="Corrente Fonte", color="grey")
plt.xlabel("Tempo de Simulacao (seg)", fontsize=12)
plt.ylabel("Amplitude [A]", fontsize=12)
plt.legend(frameon=True, facecolor="white", edgecolor="white")
plt.grid()


plt.show()



