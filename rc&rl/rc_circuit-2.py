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

Simulacao de circuito eletrico modelando resistor e indutor em serie (Circuito RC serie)

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
C = 1.0e-6  # capacitancia (Farads)

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

# Fonte de tensao
vsource = 1 * np.sin(w*t)

for time in range(0, t.shape[0]):
    if (t[time] < 0.005):
        vsource[time] = 0

########################################################################################################################
### ... ::: VALORES INICIAIS ::: ... ###
## Valor inicial das variaveis de estado (vC)
iC_ini = 0  # chute inicial: corrente do capacitor t(0-)

out_iC = [iC_ini]  # lista que armazena e atualiza os valores de iC


ks = 0.05
Rs = ks * (deltaT / (2 * C))

### ... ::: FIM VALORES INICIAIS ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: INTEGRACAO NUMERICA ::: ... ###
# for loop iniciando do ponto seguinte ao valor inicial
for n in range(1, nptos):

    if (Rs == 0):
        iC = (2 * C / deltaT) * (vsource[n] - vsource[n-1]) - iC_ini

    else:
        iC = (2 * C / (deltaT + 2 * Rs * C)) * (vsource[n] - vsource[n-1] + ((2 * Rs * C - deltaT)/ (2 * C)) * iC_ini)

    iC_ini = iC
    out_iC.append(iC_ini)

out_iC = np.array(out_iC)

### ... ::: FIM INTEGRACAO NUMERICA ::: ... ###
########################################################################################################################


########################################################################################################################
### ... ::: PLOT DOS RESULTADOS ::: ... ###
pwd = os.getcwd()

fig, ax = plt.subplots(2, 1, figsize=(9, 6))

plt.axes(ax[0])
plt.title('Resposta da Corrente no Capacitor \n para Fonte de Tensao de Excitacao Senoidal \n ks = {}'.format(ks))
plt.plot(1000 * t, vsource, '-b', label='Fonte de Tensao [V]')
plt.ylabel('Amplitude [V]')
plt.grid()
plt.legend(frameon=True, edgecolor='blue', framealpha=1.)

plt.axes(ax[1])
plt.plot(1000 * t, out_iC, '-g', label='Corrente no Capacitor [A]')
plt.ylabel('Amplitude [A]')
plt.xlabel('Tempo [ms]')
plt.grid()
plt.legend(frameon=True, edgecolor='green', framealpha=1.)

plt.savefig(pwd + '/Resposta_CircuitoRCserie_inTensao_outCorrente_ks.png', format='png')

### ... ::: FIM PLOT DOS RESULTADOS ::: ... ###
########################################################################################################################
