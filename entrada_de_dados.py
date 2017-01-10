# -*- coding: utf-8 -*-
"""
UNIVERSIDADE FEDERAL DE PERNAMBUCO
Autor: Bruno Sampaio Alves
Disciplina: Metodo dos Elementos Finitos
Professor: Paulo Marcelo Vieira Ribeiro

            Entrada de dados
"""
from elementos import ElementoCST, Barra, ElementoQ4, ElementoQ8, ElementoLST
import numpy as np  

def importar_pontos(local_arquivo):
    
    # Funcao para importar as coordenadas dos pontos a partir de um arquivo txt                                                    
    
    ##############  O ARQUIVO IMPORTADO SO PODE CONTER NUMEROS ################
    a = np.loadtxt(local_arquivo, 
                   dtype={'names': ('pt', 'x', 'y', 'z'),
                   'formats': ('i3', 'f11', 'f11','f11')})

    pontos = []               
    for i in range(len(a)):
        
        #Adiciona-se os pontos carregados a lista de pontos
        pontos.append([a[i][1],a[i][2],a[i][3]])

    return pontos

def importar_barras(local_arquivo,A,E,Iz,Iy,G,J,pontos,p):
    u"""Importa as conectividades das barras de um arquivo txt   

    O arquivo deve esta escrito no seguinte formato:
    Nº do elemento , Ponto Inicial, Ponto final         

    Note aqui que os pontos inicial e final devem ter coordenada comecando de
    uma unidade e nao de zero como normalmente é em python.                                       
    """
    ##############  O ARQUIVO IMPORTADO SO PODE CONTER NUMEROS ################
    barra1 = np.loadtxt(local_arquivo, 
                   dtype={'names': ('elemento', 'P1', 'P2'),
                   'formats': ('i3', 'i3', 'i3')})    
    
    barra = []
    
    for i in range(len(barra1)): 
        
        barra.append(Barra(pontos[barra1[i][1]-1], 
                           pontos[barra1[i][2]-1], A, E, Iz, Iy, G, J,p))

    return barra
                                   
def importar_elementos_cst(local_arquivo,E,ni,t,pontos):
    
    #Importa as conectividades das barras de um arquivo txt                                             
    
    ##############  O ARQUIVO IMPORTADO SO PODE CONTER NUMEROS ################
    CST1 = np.loadtxt(local_arquivo, 
                   dtype={'names': ('elemento', 'P1', 'P2', 'P3'),
                   'formats': ('i3', 'i3', 'i3','i3')})    
    
    CST = []
    
    for i in range(len(CST1)): 
        
        CST.append(ElementoCST(pontos[CST1[i][1]-1], 
                                pontos[CST1[i][2]-1], 
                                pontos[CST1[i][3]-1], E, ni, t))
                                   
    return CST

def importar_elementos_Q4(local_arquivo,E,ni,t,pontos):
    
    #Importa as conectividades das barras de um arquivo txt                                         
    
    ##############  O ARQUIVO IMPORTADO SO PODE CONTER NUMEROS ################
    CST1 = np.loadtxt(local_arquivo, 
                   dtype={'names': ('elemento', 'P1', 'P2', 'P3','P4'),
                   'formats': ('i3', 'i3', 'i3','i3','i3')})    
    
    CST = []
    
    for i in range(len(CST1)): 
        
        CST.append(ElementoQ4(pontos[CST1[i][1]-1], 
                                      pontos[CST1[i][2]-1], 
                                      pontos[CST1[i][3]-1], 
                                      pontos[CST1[i][4]-1], E, ni, t))
                                   
    return CST

def importar_elementos_Q8(local_arquivo,E,ni,t,pontos):
    
    #Importa as conectividades das barras de um arquivo txt                                               
    
    ##############  O ARQUIVO IMPORTADO SO PODE CONTER NUMEROS ################
    CST1 = np.loadtxt(local_arquivo, 
                   dtype={'names': ('elemento', 'P1', 'P2', 'P3','P4','P5','P6','P7','P8'),
                   'formats': ('i3', 'i3', 'i3','i3','i3','i3', 'i3','i3','i3')})    
    
    CST = []
    
    for i in range(len(CST1)): 
        
        CST.append(ElementoQ8(pontos[CST1[i][1]-1], 
                                      pontos[CST1[i][2]-1], 
                                      pontos[CST1[i][3]-1], 
                                      pontos[CST1[i][4]-1], 
                                      pontos[CST1[i][5]-1], 
                                      pontos[CST1[i][6]-1], 
                                      pontos[CST1[i][7]-1], 
                                      pontos[CST1[i][8]-1], E, ni, t))
                                   
    return CST   

def importar_elementos_LST(local_arquivo,E,ni,t,pontos):
    
    #Importa as conectividades das barras de um arquivo txt
    
                                                        
    
    ##############  O ARQUIVO IMPORTADO SO PODE CONTER NUMEROS ################
    CST1 = np.loadtxt(local_arquivo, 
                   dtype={'names': ('elemento', 'P1', 'P2', 'P3','P4','P5','P6'),
                   'formats': ('i3', 'i3', 'i3','i3','i3','i3', 'i3')})    
    
    CST = []
    
    for i in range(len(CST1)): 
        
        CST.append(ElementoLST(pontos[CST1[i][1]-1], 
                               pontos[CST1[i][2]-1], 
                               pontos[CST1[i][3]-1], 
                               pontos[CST1[i][4]-1], 
                               pontos[CST1[i][5]-1], 
                               pontos[CST1[i][6]-1], E, ni, t))
                                   
    return CST                                