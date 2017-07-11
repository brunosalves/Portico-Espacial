#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 21:17:45 2017

@author: bruno
"""
from __future__ import division
import numpy as np
from processamento import definir_fg, calcular_kg, atribuir_restricoes,\
calcular_mg, resolver_sistema, expandir_xgr
from sympy import Symbol, lambdify, Matrix
from settings import ngl, elem_matriz

class Sistema(object):
    """Classe responsavel por agrupar todas as informacoes do sistema 
    estrutural em questao.
    
    Metodos:
        adicionar_forcas
    """

    def __init__(self, pontos, barras=[], restricoes=[]):

        self.pontos = pontos
        self.n_pontos = np.shape(pontos)[0]
        self.barras = barras
        self.rest = restricoes
        self.carregamentos = []
        self.resultados = []

    def adicionar_forcas(self,lista_forcas):
        """Transforma a lista que contem as informacoes das forcas no formato
        a seguir:
            
                [ponto de aplicacao, grau de liberdade associado, valor da forca]
            
            Em um vetor de forcas para a analise por meio do metodo dos 
        elementos finitos.
        
        Input:
            lista_forcas = list
                Lista contendo as forcas nodais aplicadas 
        """
        assert np.shape(lista_forcas)[1] == 3, "Valor de entrada no formato incorreto!"
        fg =  definir_fg(lista_forcas, self.n_pontos)
        self.carregamentos.append(fg)
        self.carregamentos_r = [
                atribuir_restricoes(fg, self.rest) for fg in self.carregamentos
                ]

    def definir_mg(self):
        """Funcao que define a matriz de massa do sistema
        """
        return calcular_mg(self.pontos, self.barras)
        
    def matrizes_de_rigidez(self):
        """Funcao que define as matrizes de rigidez do sistema
        """
        self.kg = calcular_kg(self.pontos, self.barras)

        self.kgr = atribuir_restricoes(self.kg, self.rest)

    def reduzir_vetores_carregamento(self, lista_forcas):

        self.carregamentos_r = [
                atribuir_restricoes(fg, self.rest) for fg in self.carregamentos
                ]

    def analise_estatica(self):

        assert elem_matriz in [np.float32, np.float64], "A variavel elem_matriz no arquivo settings deveria esta definida como float!"
        self.matrizes_de_rigidez()
        
        xgr = (resolver_sistema(self.kgr, fgr) for fgr in self.carregamentos_r)

        x_g = (expandir_xgr(xgri, self.pontos, self.rest) for xgri in xgr)

        self.resultados = [Resultados(xgi) for xgi in x_g]
        
        for result in self.resultados:
            result.f_i = forcas_internas(self.barras, result.x_g)
#        assert all(len(result.f_i(self.barras) for result in self.resultados)), "O numero de resultados nao pode ser diferente do numero de barras!"

#%%
def forcas_internas(barras, x_g, var_i = np.nan):
    """Calcula o vetor de forcas internas referentes a cada elemento da lista
    de barras"""
    assert len(x_g) > 1, "Essa formulacao desenvolvida para analisar os vetores X cada um por vez."

    # Indices do vetor x_i associados com os graus de liberdade de cada barra
    ref = ([l for k in barra.ref_pontos for l in range(k,k+ngl)] \
            for barra in barras)
    
    # Deslocamentos associados aos nos das barras
    x_i = (x_g[i] for i in ref)

    # Forcas internas em cada barra
    if isinstance(var_i, list) or isinstance(var_i, np.ndarray):
        f_ig = [np.dot(barra.Kl(*var_i), xi) for barra in barras for xi in x_i]
    else:
        f_ig = [np.dot(barra.Kl, xi) for barra in barras for xi in x_i]

    f_i = [np.dot(barras[i].T, f_ig[i]) for i in range(len(barras))]

    assert len(barras) == np.shape(f_i)[0], "O numero de esforcos obtidos deve ser igual ao numero de barras!"

    return f_i
#%%
class SistemaOtimizacao(Sistema):
    """ Classe que representa o sistema estrutural a ser otimizado. Ate o
    momento apenas foi implementado o caso onde as variaveis de projeto sao
    as dimensoes da peca e a area de aco.
    """
    def __init__(self, variaveis, pontos, barras=[], restricoes=[]):

        assert isinstance(elem_matriz, object), "A variavel elem_matriz no arquivo settings deveria esta definida como object!"
        assert all(x.__class__ == Symbol for x in variaveis),\
            "As variaveis devem ser simbolicas!"
        
        self.var = variaveis
        super(SistemaOtimizacao, self).__init__(pontos, barras=barras,\
             restricoes=restricoes)

        self.var_i = np.array([])

        kg = calcular_kg(self.pontos, self.barras)
        kgr = atribuir_restricoes(kg, self.rest)
        self.kgr = lambdify(self.var, Matrix(kgr), "numpy")
        self.kg = lambdify(self.var, Matrix(kg), "numpy")
        
        for barra in self.barras:
            barra.Kl = lambdify(self.var, Matrix(barra.Kl), "numpy")

    def adicionar_forcas(self,lista_forcas):
        
        super(SistemaOtimizacao, self).adicionar_forcas(lista_forcas)
        
        self.carregamentos = [fg.astype(float) for fg in self.carregamentos]
        self.carregamentos_r = [
                atribuir_restricoes(fg, self.rest) for fg in self.carregamentos
                ]

    def resolver_sistema(self, var_i, calc_fi = True):
        """Resolve o sistema para um dado valor de var_i. Onde var_i trata-se
        de um vetor contendo todas as variaveis de projeto.
        """
        assert not np.array_equal(self.var_i, var_i), "Os resultados ja foram calculados!"
        self.var_i = var_i

        assert len(self.carregamentos_r) > 0, "os carregaentos devem esta definidos!"
        assert all([isinstance(fg, np.ndarray) for fg in self.carregamentos]),"fg do tipo errado!"

        kgr = self.kgr(*self.var_i)
#        x_gr_i = (resolver_sistema(kgr, fgr) for fgr in self.carregamentos_r)
        x_gr_i = (resolver_sistema(kgr, x_ri) for x_ri in self.carregamentos_r)
        _xg_i = (expandir_xgr(xgri, self.pontos, self.rest).astype(float) for
                 xgri in x_gr_i)
        
        self.resultados = [Resultados(xgi) for xgi in _xg_i]
        
        if calc_fi == True:
            for result in self.resultados:
                result.f_i = forcas_internas(self.barras, result.x_g, var_i=self.var_i)
            assert all(len(result.f_i) == len(self.barras) for result in self.resultados), "O numero de resultados nao pode ser diferente do numero de barras!"

    def obter_x(self, var_i, calc_fi= True):
        assert len(self.resultados) == 1, "A lista de resultados deveria conter apenas 1 elemento!"
        if not np.array_equal(self.var_i, var_i):
            self.resolver_sistema(var_i, calc_fi)
        return self.resultados[0].x_g
        
    def obter_dz(self, var_i, calc_fi= True):
        x = self.obter_x(var_i, calc_fi)
        return x[2:-1:6]
    
    def obter_fi1(self, var_i, calc_fi= True):
        """Retorna os esforcos axiais atuantes em cada barra.
        """
        assert len(self.resultados) == 1, "A lista de resultados deveria conter apenas 1 elemento!"
        if not np.array_equal(self.var_i, var_i):
            self.resolver_sistema(var_i, calc_fi)
        return [[f_i[0][0], f_i[6][0]] for f_i in self.resultados[0].f_i]

class Resultados(object):
    
    def __init__(self, x_g):
        assert isinstance(x_g, list) or isinstance(x_g, np.ndarray), "x_g nao esta definido com o tipo correto!"
        self.x_g = x_g # Deslocamentos dos nos do sistema
        self.f_i = [] # Lista com as forcas internas de cada elemento