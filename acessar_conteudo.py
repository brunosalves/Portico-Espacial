# -*- coding: utf-8 -*-
"""
UNIVERSIDADE FEDERAL DE PERNAMBUCO
Autor: Bruno Sampaio Alves
Disciplina: Metodo dos Elementos Finitos
Professor: Paulo Marcelo Vieira Ribeiro

################## funcoes para o usuario acessar os resultados ###############
"""

from settings import ngl
from elementos import Barra, ElementoCST

def inspecionar_ponto(x,fg):
    
    ponto = input("Qual ponto deseja-se inspecionar os resultados?\n")
    
    while (ponto > len(x)/ngl) or (ponto < 1):
        
        print "\nERRO! Este ponto nao existe. Tente novamente."
        ponto = input("De qual ponto deseja-se inspecionar os resultados?\n")
    
    print "\nDeslocamentos:"
    
    for i in range(ngl):
        
        print "x%1.0i = %s" % (i+1,x[(ponto-1)*ngl+i])
        
    print "\nForcas:"
    
    for i in range(ngl):
        
        print "F%1.0i = %s" % (i+1,fg[(ponto-1)*ngl+i])

def inspecionar_elemento(lista_elementos):

    elemento = input("Qual elemento deseja-se inspecionar os esforcos?\n")
    
    while (elemento > len(lista_elementos)) or (elemento < 1):
        
        print "\nERRO! Este elemento nao existe. Tente novamente."
        elemento = input("Qual elemento deseja-se inspecionar os esforcos?\n")
        
    print "\nEsforcos internos no elemento %s:" % (elemento)
    
    for i in range(len(lista_elementos[elemento-1].Fl)):
        
        print "f%1.0i = %s" % (i+1, lista_elementos[elemento-1].Fl[i])
    

def inspecionar(x,fg,lista_elementos,malha):
    
    print "\nResultados obtidos com sucesso!"
    
    escolha = '1'
    
    while escolha == '1' or escolha == '2':
        
        print "\nQuais resultados voce deseja consultar?"
        print "1 - Deslocamentos, forcas e reacoes nos pontos."
        print "2 - Esforcos internos em um dado elemento."
        
        escolha = raw_input("Escolha uma das opcoes acima.\n")
        
        if escolha == '1':
            
            inspecionar_ponto(x,fg)
            
        elif escolha == '2':
            
            if malha < 5:
                inspecionar_elemento(lista_elementos)
            else:
                print "Opcao nao disponivel para a malha escolhida."
        else:
            
            pass


if __name__ == "__main__":

    #inspecionar_ponto([1,2,3,4,5,6,7,8,9,10],[5,6,7,8,9,10,11,12,13,14])
    inspecionar(1,1,1)