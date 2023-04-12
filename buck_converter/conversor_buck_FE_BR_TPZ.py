# !/usr/bin/env python3
# -*- coding_ctrl: utf-8 -*-

#########################################
# Created by: Joao Pedro Peters Barbosa #
#                                       #
# email: joao.peters@engenharia.ufjf.br #
#########################################


"""
Disciplina [210081] - Tec Simulacao de Conversores Estaticos

Desenvolvimento do programa teste apresentado na aula de 06 de Maio de 2021

Simulacao de circuito eletrico modelando conversor buck

Metodos de Integracao utilizados: FE, BE, Trapezoidal
"""

########################################################################################################################
### ... ::: BIBLIOTECAS UTILIZADAS ::: ... ###
import matplotlib.pyplot as plt
import numpy as np
import os

### ... ::: FIM BIBLIOTECAS UTILIZADAS ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: DADOS DO CIRCUITO ::: ... ###
Vd = 8  # valor de tensao da fonte (volts)
rL = 1.0e-3  # resistencia serie (ohms)
L  = 5.0e-6  # indutancia serie (Henry)
C  = 100.0e-6  # capacitancia paralela (Faraday)
R = 1.0  # ohms

f = 60  # frequencia eletrica (Hertz)
w = 2 * np.pi * f  # frequencia eletrica angular (rad/s)

### ... ::: FIM DADOS DO CIRCUITO ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: DADOS DE INTEGRACAO NUMERICA ::: ... ###
## Tempo de simulacao
t_ini = 0  # segundos
t_end = 1.0e-3  # segundos

## Passo de simulacao
deltaT = 0.2e-6  # segundos

## Vetor tempo de simulacao
t = np.arange(t_ini, t_end + deltaT, deltaT)
nptos = len(t)

### ... ::: FIM DADOS DE INTEGRACAO NUMERICA ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: DADOS CONVERSOR BUCK ::: ... ###
fs = 100e3  # frequencia de chaveamento do conversor buck (Hertz)
ws = 2 * np.pi * fs  # frequencia angular de chaveamento do conversor buck (rad/s)
Ts = 1 / fs  # periodo de chaveamento do conversor buck (segundo)

## Valores de Razao Ciclica
D = 0.8  # relacao ton/t_total
Vi_ini = 0  # chute inicial valor de tensao na saida

Vi = Vi_ini

out_D = [D]  # lista que armazena valores de D
out_Vi = [Vi_ini]  # lista que armazena valores de tensao na saida

### ... ::: FIM DADOS CONVERSOR BUCK
########################################################################################################################


########################################################################################################################
### ... ::: VALORES INICIAIS ::: ... ###
## Valor inicial das variaveis de estado (iL, vC)
x_ini = np.array([0, 0])  # chute inicial

x = x_ini

out_iL = [x_ini[0]]  # lista que armazena e atualiza os valores de iL
out_vC = [x_ini[1]]  # lista que armazena e atualiza os valores de vC

## Valor inicial da onda dente de serra
saw_ini = 0
out_saw = [saw_ini]  # lista que armazena valores da onda dente de serra (comparador para funcionamento correto do duty)

### ... ::: FIM VALORES INICIAIS ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: INTEGRACAO NUMERICA ::: ... ###
# for loop iniciando do ponto seguinte ao valor inicial
for n in range(1, nptos):

    # onda dente de serra de frequencia fs
    saw = saw_ini + (deltaT / Ts)

    if saw > 1:
        saw_ini = saw - 1

    else:
        saw_ini = saw

    # tensao chaveada do conversor buck
    if D >= saw:
        Vi = Vd;

    else:
        Vi = 0

    Vi_ini = Vi

    # matrizes de estado do conversor buck
    A = np.array([[-rL / L, -1 / L], [1 / C, -1 / (R * C)]])
    B = np.array([1 / L, 0])
    I2 = np.eye(2)

    # novas matrizes de estado para modo de conducao descontinuo do indutor
    if x[0] <= 0:
        A = np.array([[0, 0], [0, -1 / (R * C)]])

    #
    # Meotodo Forward Euler
    M = I2 + (A * deltaT)
    x = np.matmul(M, x_ini) + (B * deltaT * Vi)
    flag = 1

    # #
    # # Metodo Backward Euler
    # M = np.linalg.inv(I2 - (A * deltaT))
    # x = np.matmul(M, x_ini) + np.matmul(M, B * deltaT) * Vi
    # flag = 2

    # #
    # # Metodo Trapezoidal
    # M = np.matmul(np.linalg.inv(I2 - (A * (deltaT / 2))), I2 + (A * (deltaT / 2)))
    # N = np.matmul(np.linalg.inv(I2 - (A * (deltaT / 2))), B * deltaT)
    # x = np.matmul(M, x_ini) + N * ((Vi + Vi_ini) / 2)
    # flag = 3

    # atualizacao dos valores iniciais
    x_ini = x

    # saida dos dados
    out_saw.append(saw)
    out_D.append(D)
    out_iL.append(x[0])
    out_vC.append(x[1])
    out_Vi.append(Vi)

out_saw = np.array(out_saw)
out_D = np.array(out_D)
out_iL = np.array(out_iL)
out_vC = np.array(out_vC)
out_Vi = np.array(out_Vi)

### ... ::: FIM INTEGRACAO NUMERICA ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: PLOT DOS RESULTADOS ::: ... ###
pwd = os.getcwd()

fig, ax = plt.subplots(1, figsize=(9, 6))
plt.plot(1000 * t, out_iL, 1000 * t, out_vC)

# for spine in plt.gca().spines.values():
#     spine.set_visible(False)

if flag == 1:
    metodo = 'Metodo de Integracao Numerica Forward-Euler'

elif flag == 2:
    metodo = 'Metodo de Integracao Numerica Backward-Euler'

elif flag == 3:
    metodo = 'Metodo de Integracao Numerica Trapezoidal'

plt.title('Solucao aproximada para Conversor Buck \n ' + metodo)
plt.legend(['iL', 'vC'])
plt.xlabel('Tempo [ns]')
plt.ylabel('[V] & [A]')
plt.grid()

if flag == 1:
    plt.savefig(pwd + '/ConvBuck_FE.png', format='png')

elif flag == 2:
    plt.savefig(pwd + '/ConvBuck_BE.png', format='png')

elif flag == 3:
    plt.savefig(pwd + '/ConvBuck_TPZ.png', format='png')
### ... ::: FIM PLOT DOS RESULTADOS ::: ... ###
########################################################################################################################
