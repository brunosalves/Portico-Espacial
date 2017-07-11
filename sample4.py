# -*- coding: utf-8 -*-
"""
Created on Sat Feb 04 14:35:29 2017

@author: Bruno1

Exemplo: Viga engastada com 10m de comprimento
"""
from __future__ import division
from elementos import Barra
from processamento import calcular_kg, atribuir_restricoes, definir_fg,\
analise_estatica, gerar_matrizes
from posprocessamento import mostrar_barras, inspecionar_ponto
import numpy as np
import sympy as sp
from settings import ngl

if __name__ == "__main__":
    
    # Entrada de dados
    E = 2.0 * 10**11 # Modulo de Elasticidade em Pa
    ni = 0.2 # coeficiente de poisson
    P = 10**5 # Forca atuante em Newtons
    A = 0.5 ** 2 # Area da secao transversa
    G = E/(2.0*(1.0 + ni)) # Modulo de elasticidade transversal em Pa
    J = 2.25 * ( 0.5 / 2.0 ) ** 4.0 # Constante torsional para secao retangular
    p = 0.0*A # densidade por metro linear do elemento
    I = 0.5 ** 4.0 / 12.0 # Momentos de inercia em m^4
    m_concentrada = 0
    
    # Introducao dos pontos
    pontos = [[0,0,0],[4,0,0],[4,4,0],[4,4,4]]
    
    # Definicao dos graus de liberdades restringidos na matriz de rigidez
    restricoes = [0,1,2,3,4,5]
    
    # Introduz-se as barras de cada pavimento
    barras = []
    
    for i in range(len(pontos)-1):
        #print i, i+1, pontos[i], pontos[i+1]
        barras.append( Barra( pontos[i], pontos[i+1],A,E,I,I,G,J,p))
    
    # Aplicacao das cargas
    forcas_estaticas = [ [ len(pontos), 3, -P ] ]
    
    # Passo 1 - Calcular as matrizes do sistema
    kg,kgr,mg,mgr,fgr,fg = gerar_matrizes(pontos,barras,restricoes,forcas_estaticas)  
    
    # Passo 2 - Analise estatica
    x,fg_c_reacoes = analise_estatica(kg,kgr,fgr,pontos, restricoes,forcas_estaticas)
    
    mostrar_barras(barras)
    
    inspecionar_ponto( 3, pontos, x, fg_c_reacoes )