# -*- coding: utf-8 -*-
"""
UNIVERSIDADE FEDERAL DE PERNAMBUCO
Autor: Bruno Sampaio Alves
Professor: Paulo Marcelo Vieira Ribeiro
    
        PROGRAMA PARA CALCULO DE ESFORCOS EM PORTICOS ESPACIAIS

NOMENCLATURA:    

pontos = lista com todos os pontos
barras = lista com todas as barras
elementos_cst = lista com todos os elementos cst    
forcas = lista com as forcas nodais aplicadas na estrutura. Segue o seguinte
         formato :
         [ponto de aplicacao > 0, direcao de aplicacao, valor da forca]    
cargas_distribuidas = lista com os carregamentos distribuidos aplicados na
         estrutura. Segue o seguinte formato:
         [barra submetida ao carregamento > 0, direcao de aplicacao, valor do carregamento]

kg = matriz de rigidez global do sistema
kgr = matriz de rigidez global com as linhas e colunas referentes as
      restricoes removidas
mg = matriz de massa global do sistema
mgr = matriz de massa global com as linhas e colunas referentes as
      restricoes removidas
fg = vetor de forcas do sistema
fgr = vetor de forcas do sistema com as reacoes de apoio removidas
xgr = vetor de deslocamentos obtido a partir de kgr e fgr
x = vetor de deslocamentos global incluindo os deslocamentos dos apoios (0)
fg_c_reacoes = vetor de forcas do sistema agora incluindo as reacoes de apoios
w = frequencias naturais da estrutura
v = matriz com os autovetores resultantes da analise dinamica
rest = vetor com os indices da matriz que serao restritos (indice zero \
       refere-se a linha 1 da matriz)

FORMATO DOS DADOS A SEREM ENTRADOS:

    pontos = [ [0,0,0], [1,0,0], [0,1,0], ....]
    conectividade = [[coordenada_ponto1,coordenada_ponto2], ... ]
    rest = [0,1,2, 5,6,7,...] -> cada numero refere-se a uma linha/colona
                                 da matriz de rigidez e massa global
    forcas = [[ponto_de_atuacao>0,grau_de_liberdade,valor_da_forca], ... ]

    Convencao Graus de Liberdade:
    1 = deslocamento na direcao X
    2 = deslocamento na direcao Y
    3 = deslocamento na direcao Z
    4 = rotacao em torno de X
    5 = rotacao em torno de Y
    6 = rotacao em torno de Z
"""

from settings import ngl
from processamento import calcular_kg, atribuir_restricoes, definir_fg,\
resolver_sistema, expandir_xgr, recriar_fg, calcular_esforcos_CST, \
encontrar_frequencias, calcular_mg, f_dinamico_func, \
metodo_da_diferenca_central, resposta_din_ponto_i, metodo_de_newmark
from elementos import Barra
from math import cos,pi
from entrada_de_dados import importar_barras, importar_pontos

"""
Classes e funcoes nao usadas neste trabalho:
from entrada_de_dados import importar_elementos_Q4, \
importar_elementos_cst, importar_elementos_Q8, importar_elementos_LST
from acessar_conteudo import inspecionar
"""

########################### INICIO DO PROGRAMA ################################

def gerar_matrizes(pontos,barras,restricoes,forcas,massas_pontuais):
    
    kg = calcular_kg(pontos, barras) 
    
    kgr = atribuir_restricoes(kg,restricoes)
    
    mg = calcular_mg(pontos, barras)    

    for i in massas_pontuais:
        
        mg[i[1],i[1]] = mg[i[1],i[1]] + i[0]

    mgr = atribuir_restricoes(mg,restricoes)
    
    fg = definir_fg(forcas,len(pontos))
    
    fgr = atribuir_restricoes(fg,restricoes)
    
    return kg,kgr,mg,mgr,fgr,fg

def analise_estatica(kg,kgr,fgr,pontos, restricoes, forcas):

    xgr = resolver_sistema(kgr,fgr)
    
    x = expandir_xgr(xgr,pontos, restricoes)
    
    fg_c_reacoes = recriar_fg(kgr,xgr,fgr,kg,x)  
    
    return x,fg_c_reacoes

def analise_dinamica(kgr,mgr,fgr,freq_carreg,t_carreg,dt,funcao,pontos, \
restricoes):

    w, v = encontrar_frequencias(kgr,mgr)
    
    forcas_dinamicas = f_dinamico_func(fgr,freq_carreg,t_carreg,dt,funcao)
    
    x_din,t = metodo_de_newmark(mgr,kgr,forcas_dinamicas,t_carreg,dt)

    x_din2=[]
    for i in range(len(x_din)):
        
        x_din2.append(expandir_xgr(x_din[i],pontos, restricoes))
    
    return w,v,x_din2,t,forcas_dinamicas

if __name__ == '__main__':
    
    #################### ENTRADA DE DADOS ####################################    
    #Apenas enquanto nao se providencia uma GUI ao programa
    
    
    E = 2.0 * 10**11 # Modulo de Elasticidade em Pa
    ni = 0.2 # coeficiente de poisson
    P = 1000.0 # Forca atuante em Newtons
    A = 0.5 ** 2 # Area da secao transversa
    G = E/(2.0*(1.0 + ni)) # Modulo de elasticidade transversal em Pa
    J = 2.25 * ( 0.5 / 2.0 ) ** 4.0 # Constante torsional para secao retangular
    p = 0.0*A # densidade por metro linear do elemento
    I = 0.5 ** 4.0 / 12.0 # Momentos de inercia em m^4
    m_concentrada = 500.0
    
    # Dados de entrada para se fazer a analise dinamica 
    freq_carreg = 2.0*pi/0.2 # Frequencia do carregamento aplicado em rad/s
    t_carreg = 1.0 # Tempo total que o carregamento sera aplicado
    dt = 0.001 # Passo de tempo da anlise dinamica
        
    # Geometria do problema
    pontos = importar_pontos("Portico_coordenadas.txt")
    rest = [0,1,2,3,4,5, 102,103,104,105,106,107, 108,109,110,111,112,113, \
            306,307,308,309,310,311]
    forcas_estaticas = [[26,1,P],[63,1,P]] #[62,1,P],[76,1,P]],
    barras = importar_barras("Portico_conectividades.txt",A,E,I,I,G,J,pontos,p)
    pontos_massas_pontuais = [5,9,13,21, 8,10,16,24, 14,15,23,27, 22,25,29,41,\
                              37,40,49,57, 42,45,50,61, 48,51,59,65, 58,60,66,71]    
    cargas_distribuidas = []
    massas_pontuais = []    
    # Cria uma lista onde cada elemento tem [valor da massa concentrada, /
    # grau de liberdade associado]
    for i in pontos_massas_pontuais:     
        for j in range(2):                 
            massas_pontuais.append([m_concentrada,(i-1)*ngl+j])

    # Barras sem rigidez
    sem_rigidez = range(1,81)    
    for i in sem_rigidez:
        P1 = barras[i-1].P1  
        P2 = barras[i-1].P2
        A = 0.0
        Iz = 0.0
        Iy = 0.0
        J = 0.0
        p = 0.0
        barras[i-1]=Barra(P1, P2, A, E, Iz, Iy, G, J, p)
    del P1,P2,A,Iz,Iy,J
    ########################### PROCESSAMENTO #################################     
    
    # Passo 1 - Calcular as matrizes do sistema
    kg,kgr,mg,mgr,fgr,fg = gerar_matrizes(pontos,barras,rest,forcas_estaticas,\
    massas_pontuais)  
    
    # Passo 2 - Analise estatica
    x,fg_c_reacoes = analise_estatica(kg,kgr,fgr,pontos, rest,forcas_estaticas)
    
    # Passo 3 - Analise Dinamica
    w,v,x_din,tempo,f_din = analise_dinamica(kgr,mgr,fgr,freq_carreg,\
    t_carreg,dt,cos,pontos, rest)
    
    resposta_1_ponto = resposta_din_ponto_i(378,x_din)    
    
    from matplotlib.pyplot import plot
    
    plot(tempo,resposta_1_ponto)
    print max(resposta_1_ponto),min(resposta_1_ponto)
     
    b=[]
    for i in range(len(resposta_1_ponto)):
        b.append(resposta_1_ponto[i][0])
        
    del x_din, resposta_1_ponto, f_din   
    