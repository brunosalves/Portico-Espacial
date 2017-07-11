# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 10:21:08 2017

@author: Bruno1

Universidade Federal de Pernambuco
Mestrado em Engenharia Civil - Estruturas
Disciplina: Otimizacao
Aluno: Bruno Sampaio Alves

                        Projeto final - ETAPA 02

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
expandir_xgr, analise_estatica
from posprocessamento import mostrar_barras, inspecionar_ponto
import numpy as np
from sympy import symbols, lambdify, Matrix, Symbol
from settings import ngl
from scipy.optimize import minimize
#%%
def func_volume(barras, var):
    """ Funcao que calcula internamente o volume total de concreto e a derivada
    desta funcao para implementacao no otimizador.
    """

    vol = 0

    for barra in barras:
        vol += barra.L * barra.A

    if isinstance(var, Symbol):
        vol_derivada = vol.diff(var)
    elif isinstance(var, list) or isinstance(var, tuple):
        vol_derivada = []
        for i in var:
            vol_derivada.append(vol.diff(i))
    else:
        raise Exception("Variavel de tipo desconhecido")

    vol_derivada = lambdify(var, vol_derivada, "math")

    vol = lambdify(var, vol, "math")

    volume = lambda x0: vol(*x0)
    dvolume = lambda x0: vol_derivada(*x0)

    return volume, dvolume

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

    reacoes = _fg[range(0, ngl*aux)]

    _fx = reacoes[range(0, len(reacoes), ngl)] ########## rever AUX aqui
    _fy = reacoes[range(1, len(reacoes), ngl)]
    _fz = reacoes[range(2, len(reacoes), ngl)]

    return _fx, _fy, _fz
#%%
if __name__ == "__main__":


    # Entrada de dados
    BV, HV, BP, HP = symbols("b_v h_v b_p h_p", real=True) # m Largura da secao
    X0 = [0.12, 0.5, 0.19, 0.19] # ponto inicial
    NI = 0.2 # coeficiente de poisson
    E = 2.1 * 10**7 # Modulo de Elasticidade em Pa
    G = E/(2.0*(1.0 + NI)) # Modulo de elasticidade transversal em Pa
    RO = 0.0 # densidade por metro linear do elemento
    M_CONCENTRADA = 0
    N_PAV = 2 # Numero de pavimentos da estrutura
    PE_DIREITO = 3 # m

    # Constante torsional para secao retangular
    j = lambda a, b: a * b**3 *(1/3 - 0.21*b/a *(1 - b**4/(12*a**4)))

    FZ = 100 # Forca atuante em Newtons
    FX = 10**2 # Forca atuante em Newtons

    # Pilar
    AP = BP * HP # Area da secao transversa
    IPy = BP * HP ** 3.0 / 12.0 # Momentos de inercia em m^4
    IPz = HP * BP ** 3.0 / 12.0 # Momentos de inercia em m^4
    JP = j(BP, HP)

    # Viga
    AV = BV * HV
    IVy = BV * HV ** 3.0 / 12.0
    IVz = HV * BV ** 3.0 / 12.0
    JV = j(BV, HV)

    x_max, y_max = 8, 6 # m
    dx, dy = 2, 2 # discretizacao da malha nas duas direcoes
    aux = 2 * (x_max//dx + y_max//dy) # numero de pontos no plano xy (z=0)
    auxx = x_max//dx
    auxy = y_max//dy

    # Limitadores
    P_CRz = lambdify((BP, HP),
                     pi**2 * E * IPz / (1.0 * PE_DIREITO)**2, "math")
                    # Carga critica de Euler
    P_CRy = lambdify((BP, HP),
                     pi**2 * E * IPy / (1.0 * PE_DIREITO)**2, "math")
                    # Carga critica de Euler
    DZ_MAX = min(x_max, y_max) / 250 # flecha maxima

    assert x_max % dx == 0 and y_max % dy == 0 and x_max > dx and y_max > dy,\
        "Rever modulacao!"

    # Introducao dos pontos
    PONTOS = []
    for i in range(N_PAV + 1):
        z = PE_DIREITO * i
        PONTOS += [[x, 0, z] for x in range(0, x_max+1, dx)] + \
                   [[x_max, y, z] for y in range(dy, y_max+1, dy)] + \
                   [[x, y_max, z] for x in range(x_max-dx, -1, -dx)] + \
                   [[0, y, z] for y in range(y_max-dy, 0, -dy)
                   ]

    # Definicao dos graus de liberdades restringidos na matriz de rigidez
    RESTRICOES = []
    for i in range(ngl * aux):
        RESTRICOES.append(i)

    # Introduz-se as barras de cada pavimento
    BARRAS = []
    for i in range(N_PAV):
        # Introduz os pilares de canto do pavimento 'i'
        for j in [0, # P1
                  x_max//dx, # P2
                  x_max//dx + y_max//dy, # P3
                  2 * x_max//dx + y_max//dy]: #P4
            index1 = i * aux + j
            index2 = index1 + aux
            BARRAS.append(Barra(PONTOS[index1], PONTOS[index2],
                                AP, E, IPz, IPy, G, JP, RO))
        # Introducao das vigas
        for j in range(0, aux - 1):
            index1 = i * aux + j + aux
            index2 = index1 + 1
            BARRAS.append(Barra(PONTOS[index1], PONTOS[index2],
                                AV, E, IVz, IVy, G, JV, RO))

        BARRAS.append(Barra(PONTOS[i*aux+2*aux-1], PONTOS[i*aux+aux],
                            AV, E, IVz, IVy, G, JV, RO))
    del index1, index2

    # Aplicacao das cargas
    FORCAS_VERTICAIS = [[i, 3, -FZ] for i in range(aux + 1, len(PONTOS)+1)]
    FORCAS_HORIZONTAIS = [[i, 2, FX] for i in
                          range(aux+1, len(PONTOS)+1, aux)] +\
                         [[i, 2, FX] for i in
                          range(aux+auxx+1, len(PONTOS)+1, aux)]
    FORCAS_ESTATICAS = FORCAS_VERTICAIS + FORCAS_HORIZONTAIS

    # Passo 1 - Calcular as matrizes do sistema
    KG, KGR, MG, MGR, FGR, FG = gerar_matrizes(
        PONTOS, BARRAS, RESTRICOES, FORCAS_ESTATICAS)
    del MG, MGR # nao precisa das matrizes de massa
    FGR = FGR.astype(float)
#%%
    # Transforma as expressoes literais em funcao numerica para acelerar
    # os calculos
    KGR = lambdify((BV, HV, BP, HP), Matrix(KGR), "numpy")
    KG = lambdify((BV, HV, BP, HP), Matrix(KG), "numpy")

    # Cria-se a funcao objetiva e a derivada da mesma
    VAR = [BV, HV, BP, HP]
    VOL, JAC = func_volume(BARRAS, VAR)

    # Limites dos valores de entrada
    bnds = ((0.12, min(x_max, y_max)/2), # base da viga
            (0.12, min(x_max, y_max)/2), # altura da viga
            (0.19, None), # base do pilar
            (0.19, None)  # Altura do pilar
           )

    # Restricoes
    CONS = (
        {'type': 'ineq',
         'fun': lambda xi: # Verifica a flecha maxima admissivel
                min(calcular_x(KGR(*xi), FGR, PONTOS, RESTRICOES)[3])
                + DZ_MAX},
        {'type': 'ineq',
         'fun': lambda xi: # Nao pode haver flambagem dos pilares
                -max(reacoes(KG(*xi), KGR(*xi), FGR, PONTOS, RESTRICOES)[2])
                + P_CRz(xi[2], xi[3])}
        )

    RESULT = minimize(fun=VOL, jac=JAC, x0=X0, method='SLSQP', bounds=bnds,
                      constraints=CONS, options={'disp':True})

    mostrar_barras(BARRAS)

    print 'x = ', RESULT.x

#%%
    X, DX, DY, DZ = calcular_x(KGR(*RESULT.x), FGR, PONTOS, RESTRICOES)

    print "\n DX             DY             DZ"
    print max(DX), max(DY), max(DZ)
    print min(DX), min(DY), min(DZ), "\n"

#%%
    # Passo 2 - Analise estatica
    
    X_SOL, FG_C_REACOES = analise_estatica(
        KG(*RESULT.x), KGR(*RESULT.x), FGR, PONTOS, RESTRICOES, FORCAS_ESTATICAS)

    inspecionar_ponto(aux*(N_PAV+1) - auxy - auxx//2-1,
                      PONTOS, X_SOL, FG_C_REACOES)
