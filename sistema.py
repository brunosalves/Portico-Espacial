# -*- coding: utf-8 -*-
"""
Created on Sun Feb 05 22:43:07 2017

@author: Bruno1

modulo destinado a definir a classe que vai representar o sistema em analise
e ira incorporar as funcoes para resolucao dos problemas nesta classe.
"""
from __future__ import division
import numpy as np

class Sistema(object):

    def __init__(self, pontos):
        
        self.pontos = pontos
        self.barras = []
        self.restricoes = []
        self.forcas_estaticas = []

if __name__ == "__main__":

    pass
