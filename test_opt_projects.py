#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  2 17:33:11 2017

@author: bruno

Script elaborado para se testar os programas de otimizacao
"""
from __future__ import division
import numpy as np
from sistema import Sistema, SistemaOtimizacao, forcas_internas
from elementos import Barra
from settings import elem_matriz, ngl
from optimization_project_02 import func_volume, calcular_x, reacoes
from math import pi
from processamento import gerar_matrizes
from sympy import symbols, lambdify, Matrix
from scipy.optimize import minimize
import pdb

def test_sistema_float():
    """Testa os resultados gerados pela classe Sistema para uma viga biapoiada
    com carga pontual aplicada no meio do seu vao.
    """
    if elem_matriz in [np.float32, np.float64]:
        
        E = 2.0 * 10**11 # Modulo de Elasticidade em Pa
        ni = 0.2 # coeficiente de poisson
        P = 1000.0 # Forca atuante em Newtons
        A = 0.5 ** 2 # Area da secao transversa
        G = E/(2.0*(1.0 + ni)) # Modulo de elasticidade transversal em Pa
        J = 2.25 * ( 0.5 / 2.0 ) ** 4.0 # Constante torsional para secao retangular
        p = 0.0*A # densidade por metro linear do elemento
        I = 0.5 ** 4.0 / 12.0 # Momentos de inercia em m^4
    
        L = 10
    
        PONTOS = [
                [0, 0, 0],
                [L/2, 0, 0],
                [L, 0, 0]
                ]
    
        BARRAS = [Barra(PONTOS[0], PONTOS[1], A, E, I, I, G, J, p),
                  Barra(PONTOS[1], PONTOS[2], A, E, I, I, G, J, p),
                ]
    
        REST = [0, 1, 2, 3, 12, 13, 14]
    
        sistema = Sistema(PONTOS, BARRAS, REST)
        
        FORCAS = [[2, 3, P]]
        
        sistema.adicionar_forcas(FORCAS)
        
        sistema.analise_estatica()
        assert np.isclose(sistema.resultados[0].x_g[8,0], P*L**3/(48*E*I)), "Deslocamento nao Ok!"
        r_apoio = sistema.resultados[0].f_i[0][2][0]
        momento = sistema.resultados[0].f_i[0][10][0]
        assert np.isclose(r_apoio, -P/2), "Reacao de apoio nao ok!"
        assert np.isclose(momento, P*L/4), "Momento fletor nao esta ok!"
    
        return sistema
#%%
def test_sistema_obj():
     
    if not elem_matriz in [np.float32, np.float64]:
        # Entrada de dados
        BV, HV, BP, HP = symbols("b_v h_v b_p h_p", real=True) # m Largura da secao
        X0 = [0.12, 0.5, 0.19, 0.19] # ponto inicial
        NI = 0.2 # coeficiente de poisson
        E = 2.1 * 10**7 # Modulo de Elasticidade em Pa
        G = E/(2.0*(1.0 + NI)) # Modulo de elasticidade transversal em Pa
        RO = 0.0 # densidade por metro linear do elemento
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
    
        # Cria-se a funcao objetiva e a derivada da mesma
        VAR = [BV, HV, BP, HP]
        VOL, JAC = func_volume(BARRAS, VAR)
    
        # Limites dos valores de entrada
        bnds = ((0.12, min(x_max, y_max)/2), # base da viga
                (0.12, min(x_max, y_max)/2), # altura da viga
                (0.19, None), # base do pilar
                (0.19, None)  # Altura do pilar
               )
    
        SISTEMA = SistemaOtimizacao(VAR, PONTOS, BARRAS, RESTRICOES)
        SISTEMA.adicionar_forcas(FORCAS_ESTATICAS)
        SISTEMA.resolver_sistema(X0)

        # Restricoes
        CONS = (
            {'type': 'ineq',
             'fun': lambda xi: # Verifica a flecha maxima admissivel
                    min(SISTEMA.obter_dz(xi))[0] + DZ_MAX},
            {'type': 'ineq',
             'fun': lambda xi: # Nao pode haver flambagem dos pilares
                    -np.max(SISTEMA.obter_fi1(xi)[0:4]) + P_CRz(xi[2], xi[3])}
            )

        RESULT = minimize(fun=VOL, jac=JAC, x0=X0, method='SLSQP', bounds=bnds,
                          constraints=CONS,
                          options={'disp':True})

        assert np.allclose(RESULT.x, [0.12, 0.5860537, 0.19, 0.43704327]), "Resultados da otimizacao invalidos!"

        print 'x = ', RESULT.x

    #%%
        X, DX, DY, DZ = calcular_x(SISTEMA.kgr(*RESULT.x), SISTEMA.carregamentos_r[0],
                                   SISTEMA.pontos, SISTEMA.rest)

        print "\n DX             DY             DZ"
        print max(DX), max(DY), max(DZ)
        print min(DX), min(DY), min(DZ), "\n"

        return SISTEMA

    else:

        return None
#%%
if __name__ == "__main__":
    
    sistema0 = test_sistema_float()
    sistema1 = test_sistema_obj()