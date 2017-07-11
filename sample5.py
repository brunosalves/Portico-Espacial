# -*- coding: utf-8 -*-
"""
Created on Sun Feb 05 21:13:15 2017

@author: Bruno1

Resolucao do exemplo 5.8 do livro First Course in the Finite Element Method
4a edicao de Daryl L. Logan (Pagina 262)

Os resultados foram todos verificados e conferem com o exemplo do livro.
"""
from __future__ import division
from elementos import Barra
from processamento import analise_estatica, gerar_matrizes
from posprocessamento import mostrar_barras, inspecionar_ponto
from settings import ngl

if __name__ == "__main__":

    # Entrada de dados
    E = 30 * 10**3 # Modulo de Elasticidade em Pa
    G = 10*10**3 # Modulo de elasticidade transversal em Pa
    J = 50 # Constante torsional para secao retangular
    I = 100 # Momentos de inercia em m^4
    A = 10 # Area da secao transversa
    RO = 0.0*A # densidade por metro linear do elemento
    
    PONTOS = [[100, 0, 0], [0, 0, 0], [100, 0, -100], [100, -100, 0],]
    
    BARRAS = []
    BARRAS.append(Barra(PONTOS[1], PONTOS[0], A, E, I, I, G, J, RO))
    BARRAS.append(Barra(PONTOS[2], PONTOS[0], A, E, I, I, G, J, RO))
    BARRAS.append(Barra(PONTOS[3], PONTOS[0], A, E, I, I, G, J, RO))

    RESTRICOES = [x for x in range(6,ngl*len(PONTOS))]

    # Aplicacao das cargas
    FORCAS_ESTATICAS = [[1, 2, -50],[1, 4, -1000]]

    # Passo 1 - Calcular as matrizes do sistema
    KG, KGR, MG, MGR, FGR, FG = gerar_matrizes(
        PONTOS, BARRAS, RESTRICOES, FORCAS_ESTATICAS)

    # Passo 2 - Analise estatica
    X, FG_C_REACOES = analise_estatica(
        KG, KGR, FGR, PONTOS, RESTRICOES, FORCAS_ESTATICAS)

    mostrar_barras(BARRAS)

    inspecionar_ponto(1, PONTOS, X, FG_C_REACOES)

    zz = BARRAS[1].Kl
