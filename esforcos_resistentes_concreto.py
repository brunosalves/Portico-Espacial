#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 11:30:34 2017

Modulo responsavel por calcular a resistencia de elementos de concreto armado

@author: bruno
"""
from __future__ import division
import numpy as np
import pdb

Es = 210000 # MPa Modulo de elasticidade do aco
gama_c = 1.4 # fator minorador da resistencia do concreto
gama_s = 1.15 # fator minorador da resistencia do aco

def mom_res_simples(A_s, b_v, fcd, fyd, E_s, altura_util):
    """Rotina para calculo do momento resistente ultimo para secoes
    retangulares com armadura simples - Ref: Milton v.01
    
    Returns:
        Mud : float
        Momento fletor ultimo resistente da secao
    """
    
    assert 5 <= fcd <= 90, "O fck deve ser dado em MPa!"
    assert 200 <= fyd < 500, "Fyd nao deve ser maior que o previsto em norma!"
    assert b_v > 5, "A largura da viga deve ser dada em cm!"
    assert altura_util > 10, "A altura util deve ser dada em cm!"
    
    fcd = fcd/10    #kN/cm2
    tcd = 0.85*fcd  #kN/cm2
    fyd = fyd/10    #kN/cm2
    E_s = E_s/10

    e_y = fyd/E_s

    # Profundidade da linha neutra na condicao balanceada
    x_b = 0.0035*altura_util/(0.0035 + e_y) 
    # Profundidade da linha neutra
    x_ = 1.25*A_s*fyd/(b_v*tcd)

    if 0 <= x_ <= x_b:
        Mud = 0.8*b_v*x_*(altura_util - 0.4*x_)*tcd/100; # kNm
    else:
        x_ = (-3.5*A_s*E_s + np.sqrt(3.5*A_s*E_s*(3.5*A_s*E_s+3200*b_v*\
                                                  altura_util*tcd)))/\
                                                  (1600*b_v*tcd)
        Mud = 0.8*b_v*x_*(altura_util-0.4*x_)*tcd/100 # kNm

    return Mud

def interacao_pilar_aprox(A_s, b_p, h_p, fyd, fcd, d_linha, N_d):
    """Formulacao aproximada de dimensionamento de pilares retangulares.
    
    A formulacao aqui adotada se refere a secao 3.4 do volume 3 do Jose Milton
    de Araujo.

    Unidades dos dados de entrada:    
        fcd => em kN/cm2
        dimensoes da secao => em cm
        esforco normal solicitante => kN
        As => cm2
        cobp => cm
        DN => mm
    """
    assert np.all(N_d >= 0), "Pilar tracionado!"
    assert 150 > b_p > 5 and 300 > h_p > 5, "As dimensoes devem esta em cm!"
    assert fcd <= 50 / gama_c, "A formulacao nao pode ser aplicada!"
    assert d_linha < h_p/2, "d_linha nao pode ser mais que a metade da altura!"
    assert 250/1.15 <= fyd <= 500/1.15, "fyd com unidade errada!"

    tcd = fcd * 0.85 # kN/m2
    v = N_d / (b_p * h_p * tcd) * 10
    w = A_s * fyd / (b_p * h_p * tcd)

    delta = d_linha / h_p

    beta = np.interp(v,\
                     [0., 0.5, 0.60, 0.70, 0.80, 0.9, 1.00],\
                     [1., 1.0, 0.93, 0.88, 0.88, 0.9, 0.93])

    if np.all(v <= 1):
        mi = (0.5 - delta) * beta * w + 0.468 * v * (1 - v)
    else:
        mi = (0.5 - delta) * beta * (w + 1 - v)

    M_res = b_p * ( h_p ** 2) * tcd * mi / 1000.0

    return M_res

def cisalhamento(V_d, Asw, b_w, altura_util, fck, fyk):
    """Funcao responsavel por calcular o esforco resistente relativo ao
    esforco cortante numa secao retangular seguindo as premissas do modelo II
    da NBR 6118:2014. Aqui foi desprezada a contribuicao favoravel da
    compressa da peca no dimensionamento do estribo.
    
    Input:
        V_d : float
            Esforco cortante de calculo
        Asw : float
            Area de aco de estribo em cm2 por metro
        b_w : float
            Largura da viga em centimetros
        h_v : float
            Altura da viga em centimetros
        fck : float (MPa)
            Resistencia caracteristica do concreto
        fyk : float (MPa)

    Output:
        Vrd2 : float (kN)
            Cortante maximo que a secao seria capaz de resistir
        Vrd3 : float (kN)
            Cortante resistido pela secao
    """
    assert b_w > 2 and altura_util > 2,\
        "As dimensoes da secao devem ser dadas em cm!"
    assert 20 <= fck <= 90, "Fck fora do previsto em norma!"

    V_d = np.abs(V_d)

    alfa = np.radians(45)
    theta = np.radians(90) 

#    fcd = fck/(gama_c * 10) # kN/cm2
    fywd = min(435.0, fyk/(gama_s * 10)) # kN/cm2

    Vrd2 = 0.54 * (1 - fck/250) * fck/gama_c * b_w * altura_util *\
        np.sin(theta)**2 * (1/np.tan(alfa) + 1/np.tan(theta)) /10

    if fck <= 50:
        fctm = 0.3 * fck ** (2/3)
    else:
        fctm = 2.12 * np.log(1 + 0.11*fck)

    Vc0 = 0.6 * b_w * altura_util * ((0.7 * (fctm)) / gama_c) / 10
    
    if V_d > Vc0:
        Vc1 = Vc0 * max((Vrd2 - V_d)/(Vrd2 - Vc0), 0)
    else:
        Vc1 = Vc0

    Vsw = Asw * 0.9 * altura_util * np.sin(alfa) * fywd * \
        (1/np.tan(alfa) + 1/np.tan(theta)) / 100

    Vrd3 = Vsw + Vc1
    print Vc0, V_d, Vc1, Vsw, Vrd3
    return Vrd2, Vrd3
