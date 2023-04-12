# !/usr/bin/env python3
# -*- coding_ctrl: utf-8 -*-

#########################################
# Created by: Joao Pedro Peters Barbosa #
#                                       #
# email: joao.peters@engenharia.ufjf.br #
#########################################


"""
Disciplina [210081] - Tec Simulacao de Conversores Estaticos

Desenvolvimento do programa para solucao de Lista enviada no dia 11 de Maio de 2021

Simulacao de circuito eletrico modelando resistor e capacitor em paralelo (Circuito RC)

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
C  = 10.0e-3  # capacitancia paralela (Faraday)
R = 1.0  # ohms

f = 60  # frequencia eletrica (Hertz)
w = 2 * np.pi * f  # frequencia eletrica angular (rad/s)

### ... ::: FIM DADOS DO CIRCUITO ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: DADOS DE INTEGRACAO NUMERICA ::: ... ###

## Constante de tempo
tau = R * C

## Tempo de simulacao
t_ini = 0  # segundos
t_end = 5 * tau  # segundos

## Passo de simulacao
deltaT = tau / 2  # segundos

## Vetor tempo de simulacao
t = np.arange(t_ini, t_end + deltaT, deltaT)
nptos = len(t)

### ... ::: FIM DADOS DE INTEGRACAO NUMERICA ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: VALORES INICIAIS ::: ... ###
## Valor inicial das variaveis de estado (iL, vC)
x_ini = 10  # chute inicial: tensao no capacitor

x = x_ini

out_vC = [x_ini]  # lista que armazena e atualiza os valores de vC

### ... ::: FIM VALORES INICIAIS ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: INTEGRACAO NUMERICA ::: ... ###
# for loop iniciando do ponto seguinte ao valor inicial
for n in range(1, nptos):

    # matrizes de estado do conversor buck
    A = np.array([-1 / (R * C)])
    # B = np.array([0])
    I2 = np.eye(1)

    # #
    # # Metodo Forward Euler
    # M = I2 + (A * deltaT)
    # x = M * x_ini
    # flag = 1

    # #
    # # Metodo Backward Euler
    # M = np.linalg.inv(I2 - (A * deltaT))
    # x = M * x_ini
    # flag = 2

    #
    # Metodo Trapezoidal
    M = np.matmul(np.linalg.inv(I2 - (A * (deltaT / 2))), I2 + (A * (deltaT / 2)))
    x = M * x_ini
    flag = 3

    # atualizacao dos valores iniciais
    x_ini = x[0][0]

    # saida dos dados
    out_vC.append(x[0][0])

out_vC = np.array(out_vC)

### ... ::: FIM INTEGRACAO NUMERICA ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: PLOT DOS RESULTADOS ::: ... ###
pwd = os.getcwd()

vC_real = 10 * np.exp(-t / R / C)
error = np.abs(out_vC - vC_real)

fig, ax = plt.subplots(1, figsize=(9, 6))

plt.figure(1, figsize=(12., 9.))
plt.plot(1000 * t, out_vC, 1000 * t, vC_real)
ax.text(t[-1] * 215, 8.65, 'Erro Médio\n{:.3f} V'.format(np.sum(error) / len(error)), color='black', fontsize='x-large', bbox=dict(facecolor='white', edgecolor='black'))

if flag == 1:
    metodo = 'Método de Integração Numérica Forward-Euler'
    acron = 'FE'

elif flag == 2:
    metodo = 'Método de Integração Numérica Backward-Euler'
    acron = 'BE'

elif flag == 3:
    metodo = 'Método de Integração Numérica Trapezoidal'
    acron = 'TPZ'

print(metodo)

plt.title('Solução aproximada para Circuito RC Paralelo \n ' + metodo)
plt.legend(['vC Solução ' + acron, 'vC Solução Exata'])
plt.xlabel('Tempo [ms]')
plt.ylabel('$v_c(t)$ [V]')
plt.grid()

plt.savefig(pwd + '/CircuitoRC_' + acron + '.png', format='png')


plt.figure(2, figsize=(12., 9.))
plt.plot(1000 * t, error)
plt.title('Erro entre Solução Recursiva e Solução Exata \n' + metodo)
plt.xlabel('Tempo [ms]')
plt.ylabel('Error [V]')
plt.grid()

plt.savefig(pwd + '/Error_CircuitoRC_' + acron + '.png', format='png')


plt.figure(3, figsize=(12., 9.))
plt.plot(1000 * t, error/vC_real * 100)
plt.title('Variação Percentual\n Solução Recursiva vs Solução Exata \n' + metodo)
plt.xlabel('Tempo [ms]')
plt.ylabel('Error [%]')
plt.grid()

plt.savefig(pwd + '/VariacaoPct_CircuitoRC_' + acron + '.png', format='png')

### ... ::: FIM PLOT DOS RESULTADOS ::: ... ###
########################################################################################################################
