# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 10:21:08 2017

@author: Bruno1
"""
from __future__ import division
from elementos import Barra
from processamento import calcular_kg, atribuir_restricoes, definir_fg,\
analise_estatica, gerar_matrizes
from posprocessamento import mostrar_barras, inspecionar_ponto
import numpy as np
import sympy as sp
from settings import ngl

#sp.init_printing()

if __name__ == "__main__":

    # Entrada de dados
    B = 0.5 # m Largura da secao
    E = 2.1 * 10**7 # Modulo de Elasticidade em Pa
    NI = 0.2 # coeficiente de poisson
    P = 10**5 # Forca atuante em Newtons
    A = B ** 2 # Area da secao transversa
    G = E/(2.0*(1.0 + NI)) # Modulo de elasticidade transversal em Pa
    J = 2.25 *(B / 2.0)** 4.0 # Constante torsional para secao retangular
    RO = 0.0*A # densidade por metro linear do elemento
    I = B ** 4.0 / 12.0 # Momentos de inercia em m^4
    M_CONCENTRADA = 0

    N_PAV = 2 # m
    PE_DIREITO = 3 # m
    
    x_max, y_max = 8,8
    dx, dy = 4, 4

    assert x_max % dx == 0 and y_max % dy == 0, "Rever modulacao!"

    # Introducao dos pontos
    PONTOS = []
    """teste = [ [x,y,0] for y in range(0,y_max+1,dy) 
        for x in range(0,x_max+1,dx) ]"""
    for i in range(N_PAV+1):
        z = PE_DIREITO * i
        PONTOS += [[0, 0, z], [4, 0, z], [8, 0, z], [8, 4, z],
                   [8, 8, z], [4, 8, z], [0, 8, z], [0, 4, z]]

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

    # Passo 2 - Analise estatica
    X, FG_C_REACOES = analise_estatica(
        KG, KGR, FGR, PONTOS, RESTRICOES, FORCAS_ESTATICAS)

    mostrar_barras(BARRAS)

    inspecionar_ponto(1, PONTOS, X, FG_C_REACOES)
