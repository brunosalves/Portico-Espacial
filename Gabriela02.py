#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 22 14:21:00 2017

@author: Bruno Sampaio Alves

Universidade Federal de Pernambuco
Mestrado em Engenharia Civil - Estruturas

Este script tem como objetivo reproduzir o exemplo 01 apresentado por Gabriella
Coelho na sua dissertacao de mestrado.
"""
from __future__ import division
from math import pi
import numpy as np
from sympy import symbols, lambdify
from scipy.optimize import minimize
from settings import ngl
from elementos import Barra
from posprocessamento import mostrar_barras, inspecionar_ponto
from sistema import SistemaOtimizacao
from esforcos_resistentes_concreto import interacao_pilar_aprox
import pdb

#%%
def custo(barras, var):
    """ Retorna o custo total da estrutura em dolares.
    Considera o volume de concreto, peso de aco e a area de forma no calculo.
    """
    BV1, HV1, BV2, HV2 = var[0], var[1], var[4], var[5]
    BP1, HP1, ASP1, BP2, HP2, ASP2 = var[8], var[9], var[10], var[11], var[12], var[13]

    _cc = 54 # USD/m3
    _cs = 0.55 # USD/kg
    _cf = 54 # USD/m2
    
    gama_s = 7850 # kgf/m3

    cost = 0
    for barra in barras:
        cost += _cc * barra.L * barra.A
        if barra.tag == "V1":
            cost += 2 * _cf * barra.L * (BV1 + HV1)
        if barra.tag == "V2":
            cost += 2 * _cf * barra.L * (BV2 + HV2)
        if barra.tag == "P1":
            cost += 2 * _cf * barra.L * (BP1 + HP1)
            cost += _cs * barra.L * ASP1 * gama_s
        if barra.tag == "P2":
            cost += 2 * _cf * barra.L * (BP2 + HP2)
            cost += _cs * barra.L * ASP2 * gama_s
#%%
    ###### FALTOU CONSIDERAR O PESO DE ACO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#%%
    cost1 = lambdify(var, cost, "math")
    cost2 = lambda x: cost1(*x)

    return cost2
#%%
def esforcos_internos(sistema, var_i, tag):
    """Calcula todas as areas de aco
    """
    assert tag in ["V1", "V2", "P1", "P2"], "tag invalida!"
    if not np.array_equal(sistema.var_i, var_i):
        sistema.resolver_sistema(var_i, calc_fi=True)

    f_i_gv1 = np.array([sistema.resultados[0].f_i[i] for i in\
                       range(len(sistema.barras)) if sistema.barras[i].tag=="V1"])
    f_i_gv2 = np.array([sistema.resultados[0].f_i[i] for i in\
                       range(len(sistema.barras)) if sistema.barras[i].tag=="V2"])
    f_i_gp1 = np.array([sistema.resultados[0].f_i[i] for i in\
                       range(len(sistema.barras)) if sistema.barras[i].tag=="P1"])
    f_i_gp2 = np.array([sistema.resultados[0].f_i[i] for i in\
                       range(len(sistema.barras)) if sistema.barras[i].tag=="P2"])

    if tag == "V1":
        return f_i_gv1
    elif tag == "V2":
        return f_i_gv2
    elif tag == "P1":
        return f_i_gp1
    else:
        return f_i_gp2

def md_mu(sistema, var_i, tag_pilar, fyd, fcd, d_linha):
    """Retorna a maior razao entre os momentos de dimensionamento e os momentos
    resistentes para um dado grupo de pilares.
    """
    assert tag_pilar in ["P1", "P2"], "Tag invalido!"
    
    if tag_pilar == "P1": #Obtem os dados para os pilares do grupo P1
        b_p = var_i[8]
        h_p = var_i[9]
        as_pilar = var_i[10]
        f_i = esforcos_internos(sistema, var_i, tag_pilar)
    else: # obtem os dados para pilares do grupo P2
        b_p = var_i[11]
        h_p = var_i[12]
        as_pilar = var_i[13]
        f_i = esforcos_internos(sistema, var_i, tag_pilar)
#    pdb.set_trace()
    N_d = f_i[:,0]
    M_d = f_i[:,4]
    fcd = fcd * 10**(-7)
    fyd = fyd * 10**(-6)
    b_p = b_p*100
    h_p = h_p*100
    M_u = interacao_pilar_aprox(as_pilar, b_p, h_p, fyd, fcd, d_linha, N_d)
    
    Md_Mu = np.max(M_d/M_u)

    return Md_Mu

#%%
def main():
#                          DEFINICAO DAS VARIAVEIS

    # Entrada de dados
    # Variaveis de projeto das vigas
    BV1, HV1, BV2, HV2, ASSUP1, ASINF1, ASSUP2, ASINF2 =\
        symbols("b_v1 h_v1 b_v2 h_v2 A_ssup1 A_sinf1 A_ssup2 A_sinf2",
                real=True)
    #Variaveis de projeto dos pilares
    BP1, HP1, BP2, HP2, AS1, AS2 = \
        symbols("b_p1 h_p1 b_p2 h_p2 A_sp1 A_sp2", real=True)
    DLINHA = 5 # cm
    NI = 0.2 # coeficiente de poisson
    E = 2.1 * 10**7 # Modulo de Elasticidade em Pa
    G = E/(2.0*(1.0 + NI)) # Modulo de elasticidade transversal em Pa
    RO = 0.0 # densidade por metro linear do elemento

    # Constante torsional para secao retangular
    j = lambda a, b: a * b**3 *(1/3 - 0.21*b/a *(1 - b**4/(12*a**4)))

    # Pilares do grupo 1
    AP1 = BP1 * HP1 # Area da secao transversa
    IPy1 = BP1 * HP1 ** 3.0 / 12.0 # Momentos de inercia em m^4
    IPz1 = HP1 * BP1 ** 3.0 / 12.0
    JP1 = j(BP1, HP1)

    # Pilares do grupo 2
    AP2 = BP2 * HP2 # Area da secao transversa
    IPy2 = BP2 * HP2 ** 3.0 / 12.0 # Momentos de inercia em m^4
    IPz2 = HP2 * BP2 ** 3.0 / 12.0
    JP2 = j(BP2, HP2)

    # Vigas do grupo 1
    AV1 = BV1 * HV1
    IVy1 = BV1 * HV1 ** 3.0 / 12.0
    IVz1 = HV1 * BV1 ** 3.0 / 12.0
    JV1 = j(BV1, HV1)

    # Vigas do grupo 2
    AV2 = BV2 * HV2
    IVy2 = BV2 * HV2 ** 3.0 / 12.0
    IVz2 = HV2 * BV2 ** 3.0 / 12.0
    JV2 = j(BV2, HV2)
#%%                         DEFINICAO DA GEOMETRIA

    # Discretizacao da malha
    DX = 4.5 # m
    DZ = 3.6 # m

    # Introducao dos pontos
    PONTOS = [[DX*i, 0, DZ*k] for k in xrange(0, 4) for i in xrange(0, 7)]

    # Introducao das barras
    BARRAS = []
    for i in xrange(0, 3): #introducao dos pilares externos
        BARRAS.append(Barra(PONTOS[i*7], PONTOS[i*7+7],
                            AP1, E, IPz1, IPy1, G, JP1, RO, tag="P1"))
        BARRAS.append(Barra(PONTOS[i*7+6], PONTOS[i*7+13],
                            AP1, E, IPz1, IPy1, G, JP1, RO, tag="P1"))
    for i in xrange(0, 3): #introducao dos pilares internos
        BARRAS.append(Barra(PONTOS[i*7+2], PONTOS[i*7+9],
                            AP2, E, IPz2, IPy2, G, JP2, RO, tag="P2"))
        BARRAS.append(Barra(PONTOS[i*7+4], PONTOS[i*7+11],
                            AP2, E, IPz2, IPy2, G, JP2, RO, tag="P2"))
    for i in xrange(1, 4):
        for j in xrange(0, 6):
            if PONTOS[7*i + j][0] >= 2*DX and PONTOS[7*i + j + 1][0] <= 4*DX:
                BARRAS.append(Barra(PONTOS[7*i + j], PONTOS[7*i + j + 1],
                                    AV2, E, IVz2, IVy2, G, JV2, RO,
                                    tag="V2"))
            else:
                BARRAS.append(Barra(PONTOS[7*i + j], PONTOS[7*i + j + 1],
                                    AV1, E, IVz1, IVy1, G, JV1, RO,
                                    tag="V1"))

    mostrar_barras(BARRAS) # Para verificar se a geometria esta Ok

#%%                 ATRIBUIR PROPRIEDADES AOS ELEMENTOS

    #Atribuir restricoes aos pontos nao usados
    REST = [x for x in xrange(42)] # Engaste nos pontos  que z = 0
    REST += [x for x in xrange(43, len(PONTOS)*ngl, 6)] # Dy = 0
    REST += [x for x in xrange(45, len(PONTOS)*ngl, 6)] # Rx = 0
    REST += [x for x in xrange(47, len(PONTOS)*ngl, 6)] # Rz = 0
    REST = np.unique(REST)

    # Carregamentos aplicados
    G = 16.5 * 10 **3 # N/m
    Q = 7.2 * 10 **3 # N/m

    # Coeficientes majoradores
    GAMA_G = 1.4
    GAMA_Q = 1.7

    Q_D = GAMA_G * G + GAMA_Q  * Q # kN/m

    # Parametros de projeto
    FCD = 23.5 * 10 ** 6 # N/m2
    FYD = 392 * 10 ** 6 # N/m2

    # Definindo-se as forcas
    F_Q = []
    for i in xrange(7, len(PONTOS)):
        F_Q.append([i, 3, -Q_D*4.5/2]) # [pt de aplicacao, grau de liberdade, valor]
        if PONTOS[i][0] == 0:
            F_Q.append([i, 5, Q_D*4.5**2/12])
        if PONTOS[i][0] == 27:
            F_Q.append([i, 5, -Q_D*4.5**2/12])

#%%                        OTIMIZACAO

    # Definicao da funcao objetivo
    VAR = [
           BV1, HV1, ASINF1, ASSUP1,   # Vigas do grupo 1
           BV2, HV2, ASINF2, ASSUP2,   # Vigas do grupo 2
           BP1, HP1, AS1,              # Pilares do grupo 1
           BP2, HP2, AS2               # Pilares do grupo 2 
    ]

    CUSTO = custo(BARRAS, VAR)

    # VALORES LIMITES PARA AS VARIÁVEIS DE PROJETO
    BVMIN = 0.20 # DIMENSAO MINIMA DA BASE DA VIGA
    BVMAX = 0.50 # DIMENSAO MAXIMA DA BASE DA VIGA
    HVMIN = 0.35 # DIMENSAO MINIMA DA ALTURA DA VIGA
    HVMAX = 0.90 # DIMENSAO MAXIMA DA ALTURA DA VIGA
    BPMIN = 0.30 # DIMENSAO MINIMA DA LARGURA DO PILAR
    BPMAX = 0.60 # DIMENSAO MAXIMA DA LARGURA DO PILAR
    HPMIN = 0.30 # ALTURA MINIMA DO PILAR
    HPMAX = 0.90 # ALTURA MAXIMA DO PILAR
    ASMINV = 2*pi*(19)**2/400 #
    ASMAXV = 6*pi*(22)**2/400 #
    ASMINP = 4*pi*(19)**2/400 #
    ASMAXP = 18*pi*(25)**2/400 #

    BNDS = (
            (BVMIN, BVMAX),   # Base das vigas 1
            (HVMIN, HVMAX),   # Altura das vigas 1
            (ASMINV, ASMAXV), # Area de aco inferior das vigas 1
            (ASMINV, ASMAXV), # Area de aco superior das vigas 1
            (BVMIN, BVMAX),   # Base das vigas 2
            (HVMIN, HVMAX),   # Altura das vigas 2
            (ASMINV, ASMAXV), # Area de aco inferior das vigas 2
            (ASMINV, ASMAXV), # Area de aco superior das vigas 2
            (BPMIN, BPMAX),   # base dos pilares 1
            (HPMIN, HPMAX),   # altura dos pilares 1
            (ASMINP, ASMAXP), # Area de aco dos pilares 1
            (BPMIN, BPMAX),   # base dos pilares 2
            (HPMIN, HPMAX),   # altura dos pilares 2
            (ASMINP, ASMAXP), # Area de aco dos pilares 2
    )

    CONS = (
            {'type': 'ineq',
             'fun': lambda(x): # Momento resistente dos pilares P1
                     md_mu(SISTEMA, x, "P1", FYD, FCD, DLINHA) - 1
            }
#            ,
#            {'type': 'ineq',
#             'fun': lambda(x):
#                     md_mu(sistema, x, "P2", FYD, FCD, DLINHA) - 1
#            }
    )

    # Pontos iniciais para a otimizacao
    X0 = [BVMIN, HVMIN, ASMINV, ASMINV]*2 + 2*[BPMIN, HPMIN, ASMINP]

#   Define-se o sistema
    SISTEMA = SistemaOtimizacao(VAR, PONTOS, BARRAS, REST)
    del VAR, PONTOS, BARRAS, REST
    SISTEMA.adicionar_forcas(F_Q)
    SISTEMA.resolver_sistema(X0)

    # Rotina de Otimizacao
    RESULTS = minimize(
                       fun=CUSTO, x0=X0, method='SLSQP',
                       bounds=BNDS,
                       constraints=CONS,
                       options={'disp':True}
                       )

    return SISTEMA

if __name__ == "__main__":

    sistema = main()
