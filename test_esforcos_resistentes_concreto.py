#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 16:13:14 2017

Testa o modulo Esforcos resistentes

@author: bruno
"""
from __future__ import division
import numpy as np
import esforcos_resistentes_concreto as resistente
from esforcos_resistentes_concreto import cisalhamento, interacao_pilar_aprox
import pytest

def test_mom_res_simples():
    resultado = resistente.mom_res_simples(3, 15, 14.3, 434.8, 200000, 36)
    assert np.isclose(resultado, 42.2923)

@pytest.mark.parametrize("V_d, Asw, b_w, altura_util, fck, fyk, experado",[
        (50, 6, 15, 36, 40, 500, (699.84, 116.6088)),
        (500, 6, 15, 36, 40, 500, (699.84, 77.4324)),
        (250, 4, 14, 40, 25, 500, (486, 67.2319)),
        (1000, 6, 15, 36, 40, 500, (699.84, 59.7659))
        ])

def test_cortante(V_d, Asw, b_w, altura_util, fck, fyk, experado):

    (vrd2, vrd3) = cisalhamento(V_d, Asw, b_w, altura_util, fck, fyk)

    assert np.isclose(vrd2, experado[0])
    assert np.isclose(vrd3, experado[1])

@pytest.mark.parametrize("A_s, b_p, h_p, fyd, fcd, d_linha, N_d, experado",[
        (6, 20, 20, 500/1.15, 30, 4, 500, 39.511),
        (6, 20, 20, 500/1.15, 30, 4, 700, 34.4366),
        (12, 25, 40, 500/1.15, 25, 5, 1500, 151.4578),
        (12, 25, 40, 500/1.15, 25, 5, 600, 158.867)
        ])
    
def test_interacao_pilar_aprox(A_s, b_p, h_p, fyd, fcd, d_linha, N_d, experado):
    resultado = interacao_pilar_aprox(A_s, b_p, h_p, fyd, fcd, d_linha, N_d)
    assert np.isclose(resultado, experado)
#def test_interacao_pilar_aprox():
#    resultado = interacao_pilar_aprox(6, 20, 20, 500/1.15, 30, 4, 500)
#    assert np.isclose(resultado, 39.511)
