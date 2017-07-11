#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 13:06:27 2017

@author: bruno
"""
from __future__ import division
from os.path import isfile
import numpy as np
from sympy import symbols
from cPickle import dump, load

def calcular_ke(E, A, L, G, J, Iz, Iy, tipo=object):
    #Constante que sera multiplicada pela matriz de rigidez adimensional
    EA_L = E * A / L
    EIz_L3 = E * Iz / ( L ** 3 )       
    EIy_L3 = E * Iy / ( L ** 3)        
    GJ_L = G * J / L
    
    #Matriz de rigidez
    Ke = np.zeros( (12,12), dtype=tipo)       
    
    #Parcela de esforco normal
    Ke[0,0] = EA_L
    Ke[6,0] = -EA_L
    Ke[0,6] = Ke[6,0]
    Ke[6,6] = Ke[0,0]
    
    #Parcela da flexao em torno de z
    Ke[1,1] = 12 * EIz_L3
    Ke[1,5] = 6 * EIz_L3 * L
    Ke[1,7] = (-1) * Ke[1,1]
    Ke[1,11] = Ke[1,5]
    
    Ke[5,1] = Ke[1,5]
    Ke[5,5] = 4 * EIz_L3 * ( L ** 2)
    Ke[5,7] = (-1) * Ke[1,5]
    Ke[5,11] = 2 * EIz_L3 * ( L ** 2 )

    Ke[7,1] = Ke[1,7]
    Ke[7,5] = Ke[5,7]
    Ke[7,7] = Ke[1,1]
    Ke[7,11] = Ke[5,7]
    
    Ke[11,1] = Ke[1,11]
    Ke[11,5] = Ke[5,11]
    Ke[11,7] = Ke[7,11]
    Ke[11,11] = 4 * EIz_L3 * ( L ** 2 )

    #Parcela da flexao em y
    Ke[2,2] = 12 * EIy_L3
    Ke[2,4] = ( -6 ) * EIy_L3 * L
    Ke[2,8] = ( -12 ) * EIy_L3
    Ke[2,10] = ( -6 ) * EIy_L3 * L
    
    Ke[4,2] = Ke[2,4]
    Ke[4,4] = 4 * EIy_L3 * ( L ** 2 )
    Ke[4,8] = 6 * EIy_L3 * L
    Ke[4,10] = 2 * EIy_L3 * ( L ** 2 )
    
    Ke[8,2] = Ke[2,8]
    Ke[8,4] = Ke[4,8]
    Ke[8,8] = Ke[2,2]
    Ke[8,10] = Ke[4,8]
    
    Ke[10,2] = Ke[2,10]
    Ke[10,4] = Ke[4,10]
    Ke[10,8] = Ke[8,10]
    Ke[10,10] = Ke[4,4]
    
    #Parcela do torsor
    Ke[3,3] = GJ_L
    Ke[3,9] = -GJ_L
    Ke[9,3] = Ke[3,9]
    Ke[9,9] = Ke[3,3]
    
    return Ke

def calcular_lambda(l, m, n, tipo=object):
    """Matriz que ira compor a matriz T
    """
    D = ( l ** 2 + m ** 2 ) ** 0.5        

    Lambda = np.zeros( ( 3 , 3 ), dtype=tipo)
    Lambda[0,0] = l
    Lambda[0,1] = m
    Lambda[0,2] = n

    if l == 0 and m == 0:
        Lambda[1,1] = 1
        Lambda[1,0] = 0
        Lambda[2,0] = (-1.0)*n
        Lambda[2,1] = 0
    else:
        Lambda[1,1] = l/D
        Lambda[1,0] = -m/D
        Lambda[2,0] = -l*n/D
        Lambda[2,1] = -m*n/D
    
    Lambda[1,2] = 0
    Lambda[2,2] = D      
    
    return Lambda

def calcular_T(Lambda, tipo=object):
    """Matriz rotacao
    """
    T = np.zeros((12, 12), dtype=tipo)

    for i in range(4):
        
        T[ i*3   , i*3   ] = Lambda[0,0]
        T[ i*3   , i*3+1 ] = Lambda[0,1]
        T[ i*3   , i*3+2 ] = Lambda[0,2]
        T[ i*3+1 , i*3   ] = Lambda[1,0]
        T[ i*3+1 , i*3+1 ] = Lambda[1,1]
        T[ i*3+1 , i*3+2 ] = Lambda[1,2]
        T[ i*3+2 , i*3   ] = Lambda[2,0]
        T[ i*3+2 , i*3+1 ] = Lambda[2,1]
        T[ i*3+2 , i*3+2 ] = Lambda[2,2]
        
    return T

def calcular_kl(E, A, L, G, J, Iz, Iy, l, m, n):
    Ke = calcular_ke(E, A, L, G, J, Iz, Iy, tipo=object)
    Lambda = calcular_lambda(l, m, n)        
    T = calcular_T(Lambda)
    # [Kl] =  [T]^t . [Ke] . [T]   
    kl = np.dot(T.T, Ke)
    kl = np.dot(kl, T) #Matriz local rotacionada para os eixos globais    
    
    return kl

def main():
    nome_arquivo = "teste.br"
    if isfile(nome_arquivo) is False:
        E, A, L, G, J, Iz, Iy, l, m, n = symbols("E A L G J Iz Iy l m n", real=True)
        kl = calcular_kl(E, A, L, G, J, Iz, Iy, l, m, n)
        dump(kl, open(nome_arquivo, "wb"))
    else:
        kl = load(open(nome_arquivo, "rb")) 
        
    return kl

if __name__ == "__main__":
    
    main()