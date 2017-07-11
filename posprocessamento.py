# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 17:29:37 2017

@author: Bruno1

Pos processamento
"""
from __future__ import division
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from settings import ngl

def mostrar_barras(lista_barras):
    """Display bars in a 3d view.
    """
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")    
    
    for barra in lista_barras:
        x,y,z = [] , [] , []
        for p_i in barra.lista_pontos:
            x_i,y_i,z_i = p_i
            x.append(x_i)
            y.append(y_i)
            z.append(z_i)    
        ax.plot_wireframe(x,y,z)
    plt.show()

def inspecionar_ponto(i,lista_pontos,x,fg):
    """Informa os esforcos pontuais existentes no ponto 'i'
    
    args (int) : numero de referencia do ponto a ser averiguado (> 0)
    """

    assert type(i) == int and i>0, "O numero de referencia do ponto deve ser \
    maior que zero."
    
    print "\nPonto %3i\n" % i
    
    print "Coordenadas ponto = %s\n" % lista_pontos[i-1]
    
    print "Dx = %s m" % x[ (i-1) * ngl ]
    print "Dy = %s m" % x[ (i-1) * ngl + 1 ]
    print "Dz = %s m" % x[ (i-1) * ngl + 2 ]
    print "Rx = %s " % x[ (i-1) * ngl + 3 ]
    print "Ry = %s " % x[ (i-1) * ngl + 4 ]
    print "Rz = %s \n" % x[ (i-1) * ngl + 5 ]
    
    print "Fx = %.3f N" % fg[ (i-1) * ngl ]
    print "Fy = %.3f N" % fg[ (i-1) * ngl + 1 ]
    print "Fz = %.3f N" % fg[ (i-1) * ngl + 2 ]
    print "Mx = %.3f Nm" % fg[ (i-1) * ngl + 3 ]
    print "My = %.3f Nm" % fg[ (i-1) * ngl + 4 ]
    print "Mz = %.3f Nm \n" % fg[ (i-1) * ngl + 5 ]
    

if __name__ == "__main__":
    pass
