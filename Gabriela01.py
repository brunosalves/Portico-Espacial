#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 16:55:14 2017

@author: Bruno Sampaio Alves

Universidade Federal de Pernambuco
Mestrado em Engenharia Civil - Estruturas

Este script tem como objetivo reproduzir o exemplo 01 apresentado por Gabriella
Coelho na sua dissertacao de mestrado.

"""
from __future__ import division
from scipy.optimize import minimize

def main():
    
    # Define-se as variaveis de projeto
    L = 5.5 #m
    alfa = 1
    w = 1.5 #tf/m
    sigma_adm = 1300 #tf/m2

    # Define-se as funcoes
    def f(x): # Funcao objetivo
        
        return 2*alfa*L*x[0]**2 + L*x[1]**2
    
    def g1(x): # Restricao 1
        
        return w*L/(2*x[0]**2) + 3*w*L**2*x[0]/(6*x[0]**4 + 4*alfa*x[1]**4) -\
                sigma_adm
    
    def g2(x): # Restricao 2
        
        return x[0]**4*(x[1] + 6*alfa*L)*w*L/\
                (x[1]**3*(6*x[0]**4 + 4*alfa*x[1]**4)*2*alfa) - sigma_adm
    
    def g3(x): # Restricao 3
        
        return (x[0]**4*x[1] + alfa*L*(3*x[0]**4 + 6*alfa*x[1]**4))*w*L/\
                (x[1]**3*(6*x[0]**4 + 4*alfa*x[1]**4)*2*alfa) - sigma_adm

    # Restricoes
    cons = ({'type': 'ineq', 'fun': lambda x: -g1(x)},
            {'type': 'ineq', 'fun': lambda x: -g2(x)},
            {'type': 'ineq', 'fun': lambda x: -g3(x)})
    
    bnds = ((0.01, 1), (0.01, 1))

    res = minimize(f, (0.5, 0.5), method='SLSQP', bounds=bnds, constraints=cons,\
                   options={'disp': True})
    
    print res.x

if __name__ == "__main__":
    
    main()