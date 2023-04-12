# !/usr/bin/env python3
# -*- coding_ctrl: utf-8 -*-

#########################################
# Created by: Joao Pedro Peters Barbosa #
#                                       #
# email: joao.peters@engenharia.ufjf.br #
#########################################


"""
Disciplina [210081] - Tec Simulacao de Conversores Estaticos

Desenvolvimento do programa para solucao de Lista2.pdf enviada no dia 20 de Maio de 2021

Simulacao de circuito eletrico modelando resistor e indutor em serie (Circuito RL serie)

Metodos de Integracao utilizados: Trapezoidal, R-K 2ª ordem, modelo de acompanhamento
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
Vs = 10  # valor de tensao da fonte (volts)

L  = 10.0e-3  # indutancia serie (Henrys)
R = 1.0  # ohms

f = 60  # frequencia eletrica (Hertz)
w = 2 * np.pi * f  # frequencia eletrica angular (rad/s)

### ... ::: FIM DADOS DO CIRCUITO ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: DADOS DE INTEGRACAO NUMERICA ::: ... ###

## Constante de tempo
tau = L / R

## Tempo de simulacao
t_ini = 0  # segundos
t_end = 10 * tau  # segundos

## Passo de simulacao
deltaT = tau / 10  # segundos

## Vetor tempo de simulacao
t = np.arange(t_ini, t_end + deltaT, deltaT)
nptos = len(t)

### ... ::: FIM DADOS DE INTEGRACAO NUMERICA ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: VALORES INICIAIS ::: ... ###
## Valor inicial das variaveis de estado (iL)
iL_ini = 0  # chute inicial: corrente do indutor t(0-)
vL_ini = 10  # chute inicial: tensao no indutor t(0-)

out_iL = [iL_ini]  # lista que armazena e atualiza os valores de iL
out_vL = [vL_ini]  # lista que armazena e atualiza os valores de vL

### ... ::: FIM VALORES INICIAIS ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: INTEGRACAO NUMERICA ::: ... ###
# for loop iniciando do ponto seguinte ao valor inicial
for n in range(1, nptos):

    # matrizes de estado do conversor buck
    A = np.array([-R / L])
    B = np.array([1 / L])
    I2 = np.eye(1)

    # #
    # # Metodo Forward Euler
    # M = I2 + (A * deltaT)
    # iL = np.matmul(M, [iL_ini]) + (B * deltaT * Vs)
    # flag = 1

    # #
    # # Metodo Backward Euler
    # M = np.linalg.inv(I2 - (A * deltaT))
    # iL = np.matmul(M, [iL_ini]) + np.matmul(M, B * deltaT) * Vs
    # flag = 2

    # #
    # # Metodo Trapezoidal
    # M = np.matmul(np.linalg.inv(I2 - (A * (deltaT / 2))), I2 + (A * (deltaT / 2)))
    # N = np.matmul(np.linalg.inv(I2 - (A * (deltaT / 2))), B * deltaT)
    # iL = np.matmul(M, [iL_ini]) + N * Vs
    # flag = 3

    # #
    # # Metodo Runge-Kutta de 2a Ordem
    # iL1 = np.matmul(I2 + (A * deltaT), [iL_ini]) + (B * deltaT * Vs)
    # diL1 = np.matmul(I2 + (A * deltaT / 2), [iL_ini]) + (B * deltaT * Vs) / 2
    # diL2 = np.matmul(A * deltaT / 2, iL1) + (B * deltaT * Vs) / 2
    #
    # iL = diL1 + diL2
    # flag = 4

    # #
    # # Metodo Runge-Kutta de 4a Ordem
    # iL1 = np.matmul(I2 + (A * deltaT) / 2, [iL_ini]) + (B * deltaT * Vs) / 2
    # iL2 = iL_ini + np.matmul(A * deltaT / 2, iL1) + (B * deltaT * Vs) / 2
    # iL3 = iL_ini + np.matmul(A * deltaT, iL2) + (B * deltaT * Vs)
    #
    # diL1 = np.matmul(I2 + (A * deltaT) / 6, [iL_ini]) + (B * deltaT * Vs) / 6
    # diL2 = np.matmul(A * deltaT / 3, iL1) + (B * deltaT * Vs) / 3
    # diL3 = np.matmul(A * deltaT / 3, iL2) + (B * deltaT * Vs) / 3
    # diL4 = np.matmul(A * deltaT / 6, iL3) + (B * deltaT * Vs) / 6
    #
    # iL = diL1 + diL2 + diL3 + diL4
    # flag = 5

    # #
    # # Modelo de Acompanhamento
    # # Aplicacao do Metodo Trapezoidal
    # i_source = vL_ini * (deltaT / 2 / L) + iL_ini
    # vL = (10 - i_source) * ((2 * L * R) / (R * deltaT + 2 * L))
    # iL = 10 - (vL / R)
    #
    # iL_ini = iL
    # vL_ini = vL
    #
    # out_iL.append(iL)
    # out_vL.append(vL)
    #
    # flag = 6

    #
    # Modelo de Acompanhamento
    # Aplicação do Metodo Runge-Kutta de 2a Ordem
    iL1 = iL_ini + (deltaT / L) * vL_ini
    vL1 = (10 - iL1) * R

    iL = iL_ini + (deltaT / 2 / L) * (vL_ini + vL1)
    vL = 10 - (iL * R)

    iL_ini = iL
    vL_ini = vL

    out_iL.append(iL)
    out_vL.append(vL)

    flag = 7




    # # atualizacao dos valores iniciais
    # iL_ini = iL[0]
    #
    # # saida dos dados
    # out_iL.append(iL[0])

out_iL = np.array(out_iL)
out_vL = np.array(out_vL)

### ... ::: FIM INTEGRACAO NUMERICA ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: PLOT DOS RESULTADOS ::: ... ###
pwd = os.getcwd()

# iL_exato = (Vs / L) * t
iL_exato = (Vs / R) * (1 - np.exp(-(R / L) * t))
error = np.abs(out_iL - iL_exato)

fig, ax = plt.subplots(1, figsize=(9, 6))

plt.figure(1, figsize=(12., 9.))
plt.plot(1000 * t, out_iL, 1000 * t, iL_exato)
ax.text(t[-1] * 215, 2.5, 'Erro Médio\n{:.3f} A'.format(np.sum(error) / len(error)), color='black', fontsize='x-large', bbox=dict(facecolor='white', edgecolor='black'))

if flag == 1:
    metodo = 'Método de Integração Numérica Forward-Euler'
    acron = 'FE'

elif flag == 2:
    metodo = 'Método de Integração Numérica Backward-Euler'
    acron = 'BE'

elif flag == 3:
    metodo = 'Método de Integração Numérica Trapezoidal'
    acron = 'TPZ'

elif flag == 4:
    metodo = 'Método de Integração Numérica Runge-Kutta de 2ª Ordem'
    acron = 'RK2'

elif flag == 5:
    metodo = 'Método de Integração Numérica Runge-Kutta de 4ª Ordem'
    acron = 'RK4'

elif flag == 6:
    metodo = 'Modelo de acompanhamento\n Aplicação de Método de Integração Numérico Trapezoidal'
    acron = 'MA-TPZ'

elif flag == 7:
    metodo = 'Modelo de acompanhamento\n Aplicação de Método de Integração Numérico Runge-Kutta de 2ª Ordem'
    acron = 'MA-RK2'

print(metodo)

plt.title('Solução aproximada para Circuito RL Serie \n ' + metodo)
plt.legend(['iL Solução ' + acron, 'iL Solução Exata'])
plt.xlabel('Tempo [ms]')
plt.ylabel('$i_L(t)$ [A]')
plt.grid()

plt.savefig(pwd + '/CircuitoRLs_' + acron + '.png', format='png')


plt.figure(2, figsize=(12., 9.))
plt.plot(1000 * t, error)
plt.title('Erro entre Solução Recursiva e Solução Exata \n' + metodo)
plt.xlabel('Tempo [ms]')
plt.ylabel('Error [A]')
plt.grid()

plt.savefig(pwd + '/Error_CircuitoRLs_' + acron + '.png', format='png')


plt.figure(3, figsize=(12., 9.))
plt.plot(1000 * t, error/iL_exato * 100)
plt.title('Variação Percentual\n Solução Recursiva vs Solução Exata \n' + metodo)
plt.xlabel('Tempo [ms]')
plt.ylabel('Error [%]')
plt.grid()

plt.savefig(pwd + '/VariacaoPct_CircuitoRLs_' + acron + '.png', format='png')

### ... ::: FIM PLOT DOS RESULTADOS ::: ... ###
########################################################################################################################
