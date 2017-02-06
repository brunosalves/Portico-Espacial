# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 13:11:37 2016
Exemplo teste de viga biapoiada
@author: Bruno1
"""
from __future__ import division
from elementos import Barra
from processamento import calcular_kg, atribuir_restricoes, definir_fg
from solucionar_problema import analise_estatica
from posprocessamento import mostrar_barras
import numpy as np

if __name__ == "__main__":

    E = 2.0 * 10**11 # Modulo de Elasticidade em Pa
    ni = 0.2 # coeficiente de poisson
    P = 1000.0 # Forca atuante em Newtons
    A = 0.5 ** 2 # Area da secao transversa
    G = E/(2.0*(1.0 + ni)) # Modulo de elasticidade transversal em Pa
    J = 2.25 * ( 0.5 / 2.0 ) ** 4.0 # Constante torsional para secao retangular
    p = 0.0*A # densidade por metro linear do elemento
    I = 0.5 ** 4.0 / 12.0 # Momentos de inercia em m^4
    m_concentrada = 500.0

    pontos = [[0,0,0],[1,0,0],[2,0,0]]
    restricoes = [2,14]
    barras = [Barra(pontos[0], pontos[1], A, E, I, I, G, J, p),
              Barra(pontos[1], pontos[2], A, E, I, I, G, J, p)
    ]
    
    forcas_estaticas = [[2,3,-1000]]
    
    kg = calcular_kg(pontos, barras) 
    
    kgr = atribuir_restricoes(kg,restricoes)
    
    fg = definir_fg(forcas_estaticas,len(pontos))
    
    fgr = atribuir_restricoes(fg,restricoes)
    
    x,fg_c_reacoes = analise_estatica(kg,kgr,fgr,pontos, restricoes,forcas_estaticas)
    
    mostrar_barras(barras)