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
L  = 1.0e-3  # indutancia serie (Henrys)

f = 50  # frequencia eletrica (Hertz)
w = 2 * np.pi * f  # frequencia eletrica angular (rad/s)

### ... ::: FIM DADOS DO CIRCUITO ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: DADOS DE INTEGRACAO NUMERICA ::: ... ###
## Tempo de simulacao
t_ini = 0  # segundos
t_end = 20e-3  # segundos

## Passo de simulacao
deltaT = 100e-6  # segundos

## Vetor tempo de simulacao
t = np.arange(t_ini, t_end + deltaT, deltaT)
nptos = len(t)

### ... ::: FIM DADOS DE INTEGRACAO NUMERICA ::: ... ###
########################################################################################################################

# Fonte de corrente
isource = 1 * np.sin(w*t)

for time in range(0, t.shape[0]):
    if (t[time] < 0.005):
        isource[time] = 0

########################################################################################################################
### ... ::: VALORES INICIAIS ::: ... ###
## Valor inicial das variaveis de estado (iL)
vL_ini = 0.0  # chute inicial: tensao no indutor t(0-)

out_vL = [vL_ini]  # lista que armazena e atualiza os valores de vL


kp = 0.
Rp = kp * (2 * L / deltaT)  # ohms

### ... ::: FIM VALORES INICIAIS ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: INTEGRACAO NUMERICA ::: ... ###
# for loop iniciando do ponto seguinte ao valor inicial
for n in range(1, nptos):

    if (Rp == 0):
        vL = (2 * L / deltaT) * (isource[n] - isource[n-1]) - vL_ini

    else:
        vL = (2 * L * Rp / (deltaT * Rp + 2 * L)) * (isource[n] - isource[n-1] + ((2 * L - Rp * deltaT)/(2 * L * Rp)) * vL_ini)

    vL_ini = vL
    out_vL.append(vL_ini)

out_vL = np.array(out_vL)

### ... ::: FIM INTEGRACAO NUMERICA ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: PLOT DOS RESULTADOS ::: ... ###
pwd = os.getcwd()

fig, ax = plt.subplots(2, 1, figsize=(9, 6))

plt.axes(ax[0])
plt.title('Resposta da Tensao no Indutor \n para Fonte de Corrente de Excitacao Senoidal \n kp = {}'.format(kp))
plt.plot(1000 * t, isource, '-g', label='Fonte de Corrente [A]')
plt.ylabel('Amplitude [A]')
plt.grid()
plt.legend(frameon=True, edgecolor='green', framealpha=1.)

plt.axes(ax[1])
plt.plot(1000 * t, out_vL, '-b', label='Tensao no Indutor [V]')
plt.ylabel('Amplitude [V]')
plt.xlabel('Tempo [ms]')
plt.grid()
plt.legend(frameon=True, edgecolor='blue', framealpha=1.)

plt.savefig(pwd + '/Resposta_CircuitoRLserie_inCorrente_outTensao_kp.png', format='png')

### ... ::: FIM PLOT DOS RESULTADOS ::: ... ###
########################################################################################################################
