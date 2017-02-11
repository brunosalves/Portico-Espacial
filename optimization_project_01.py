# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 10:21:08 2017

@author: Bruno1

Universidade Federal de Pernambuco
Mestrado em Engenharia Civil - Estruturas
Disciplina: Otimizacao
Aluno: Bruno Sampaio Alves

                        Projeto final - ETAPA 01
                        
Codigo elaborado para se otimizar a seguinte estrutura:

    Portico 3d com 3 pavimentos
    Carregamento vertical nodal nas vigas e lateral (vento) nos pilares
    neste modelo as vigas e os pilares tem a mesma secao quadrada

    minimizar:

        Volume total de concreto

    Sujeito a:
        flecha maxima <= L / 250 (NBR 6118:2014)
        reacoes verticais de apoio <= carga critica de euler
        largura minima da viga = 12 cm (NBR 6118:2014)

O objetivo desta etapa e de implementar o modulo de otimizacao na estrutura
definida na etapa 00 e estudar como ele funciona. A etapa seguinte se 
reformulara o codigo para que as secoes dos pilares e vigas sejam diferentes.

Para se executar este programa, o arquivo settings.py deve esta definido da
seguinte forma:

    ngl = 6
    elem_matriz = object

"""
from __future__ import division
from math import pi
from elementos import Barra
from processamento import gerar_matrizes, resolver_sistema,\
expandir_xgr
from posprocessamento import mostrar_barras, inspecionar_ponto
import numpy as np
from sympy import symbols, lambdify, Matrix
from settings import ngl
from scipy.optimize import minimize
#%%
def func_volume(barras, _b_):
    """ Funcao que calcula internamente o volume total de concreto e a derivada
    desta funcao para implementacao no otimizador.
    """

    vol = 0

    for barra in barras:
        vol += barra.L * barra.A

    vol_derivada = vol.diff(_b_)
    vol_derivada = lambdify(_b_, vol_derivada, "math")

    vol = lambdify(_b_, vol, "math")

    return vol, vol_derivada

def calcular_x(kgr, fgr, pontos, restricoes):
    """Resolve o sistema de equacoes formado por [kgr] * [x] = [fgr] e retorna
    o vetor de deslocamentos 'x' e os vetores contendo apenas os deslocamentos
    em cada direcao ('dx', 'dy', 'dz')
    """

    xgr = resolver_sistema(kgr, fgr)

    _x_ = expandir_xgr(xgr, pontos, restricoes).astype(float)

    _dx = _x_[range(0, len(_x_), 6)]
    _dy = _x_[range(1, len(_x_), 6)]
    _dz = _x_[range(2, len(_x_), 6)]

    return _x_, _dx, _dy, _dz

def reacoes(_kg, kgr, fgr, pontos, restricoes):
    """Calcula o vetor de forcas que incluira tambem as reacoes de apoio
    """

    _x_, _dx, _dy, _dz = calcular_x(kgr, fgr, pontos, restricoes)

    _fg = np.dot(_kg, _x_)

    reacoes = _fg[range(0, ngl*8)]

    _fx = reacoes[range(0, len(reacoes), 6)]
    _fy = reacoes[range(1, len(reacoes), 6)]
    _fz = reacoes[range(2, len(reacoes), 6)]

    return _fx, _fy, _fz
#%%
if __name__ == "__main__":


    # Entrada de dados
    B = symbols("b", real=True) # m Largura da secao
    B0 = 0.5 # ponto inicial
    E = 2.1 * 10**7 # Modulo de Elasticidade em Pa
    NI = 0.2 # coeficiente de poisson
    P = 10**5 # Forca atuante em Newtons
    A = B ** 2 # Area da secao transversa
    G = E/(2.0*(1.0 + NI)) # Modulo de elasticidade transversal em Pa
    J = 2.25 *(B / 2.0)** 4.0 # Constante torsional para secao retangular
    RO = 0.0*A # densidade por metro linear do elemento
    I = B ** 4.0 / 12.0 # Momentos de inercia em m^4
    M_CONCENTRADA = 0
    N_PAV = 3 # Numero de pavimentos da estrutura
    PE_DIREITO = 3 # m

    # Limitadores
    P_CR = lambdify(B, pi**2 * E * I / (1.0 * PE_DIREITO)**2, "math")
        # Carga critica de Euler
    DZ_MAX = N_PAV * PE_DIREITO / 250 # flecha maxima

    x_max, y_max = 8, 6
    dx, dy = 4, 3 # discretizacao da malha nas duas direcoes

    assert x_max % dx == 0 and y_max % dy == 0 and x_max > dx and y_max > dy,\
        "Rever modulacao!"

    # Introducao dos pontos
    PONTOS = []
    for i in range(N_PAV + 1):
        z = PE_DIREITO * i
        PONTOS += [[0, 0, z], [dx, 0, z], [2*dx, 0, z], [2*dx, dy, z],
                   [2*dx, 2*dy, z], [dx, 2*dy, z], [0, 2*dy, z], [0, dy, z]]

    # Definicao dos graus de liberdades restringidos na matriz de rigidez
    RESTRICOES = []
    for i in range(ngl*8):
        RESTRICOES.append(i)

    # Introduz-se as barras de cada pavimento
    BARRAS = []
    for i in range(N_PAV):
        # Introduz os pilares de canto do pavimento 'i'
        for j in [0, 2, 4, 6]:
            index1 = i*8 + j
            index2 = index1 + 8
            BARRAS.append(Barra(PONTOS[index1], PONTOS[index2],
                                A, E, I, I, G, J, RO))
        # Introducao das vigas
        for j in range(0, 7):
            index1 = i*8 + j + 8
            index2 = index1 + 1
            BARRAS.append(Barra(PONTOS[index1], PONTOS[index2],
                                A, E, I, I, G, J, RO))
        BARRAS.append(Barra(PONTOS[i*8+7+8], PONTOS[i*8+8],
                            A, E, I, I, G, J, RO))
    del index1, index2

    # Aplicacao das cargas
    FORCAS_ESTATICAS = [[10, 3, -P]]

    # Passo 1 - Calcular as matrizes do sistema
    KG, KGR, MG, MGR, FGR, FG = gerar_matrizes(
        PONTOS, BARRAS, RESTRICOES, FORCAS_ESTATICAS)
    del MG, MGR # nao precisa das matrizes de massa
    FGR = FGR.astype(float)

    # Transforma as expressoes literais em funcao numeroca para acelerar
    # os calculos
    KGR = lambdify(B, Matrix(KGR), "numpy")
    KG = lambdify(B, Matrix(KG), "numpy")

    # Cria-se a funcao objetiva e a derivada da mesma
    VOL, JAC = func_volume(BARRAS, B)

    # Restricoes
    CONS = ({'type': 'ineq', 'fun': lambda b: b - 0.12}, # b >= 12 cm
            {'type': 'ineq',
             'fun': lambda b: # Verifica a flecha maxima admissivel
                    -max(calcular_x(KGR(b[0]), FGR, PONTOS, RESTRICOES)[3])
                    + DZ_MAX},
            {'type': 'ineq',
             'fun': lambda b: # Nao pode haver flambagem dos pilares
                    -max(reacoes(KG(b[0]), KGR(b[0]),
                                 FGR, PONTOS, RESTRICOES)[2])
                    + P_CR(b[0])}
           )

    print minimize(fun=VOL, x0=B0, jac=JAC, method='SLSQP',
                   constraints=CONS, options={'disp':True})

    mostrar_barras(BARRAS)

    X, DX, DY, DZ = calcular_x(KGR(0.86568434), FGR, PONTOS, RESTRICOES)


#    # Passo 2 - Analise estatica
#    X, FG_C_REACOES = analise_estatica(
#        KG, KGR, FGR, PONTOS, RESTRICOES, FORCAS_ESTATICAS)
#
#    mostrar_barras(BARRAS)
#
#    inspecionar_ponto(1, PONTOS, X, FG_C_REACOES)
