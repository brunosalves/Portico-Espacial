# -*- coding: utf-8 -*-
"""
UNIVERSIDADE FEDERAL DE PERNAMBUCO
Autor: Bruno Sampaio Alves
Disciplina: Metodo dos elementos Finitos
Professor: Paulo Marcelo Vieira Ribeiro

                        FUNCOES DE PROCESSAMENTO
"""
from __future__ import division
import warnings
from settings import ngl, elem_matriz
from scipy.sparse import lil_matrix
#from scipy.sparse.linalg import eigsh, inv
#from scipy.linalg import solve_banded # Para solucionar um problema com sparse_matrix
import numpy as np
import threading

#%%
#Funcao para definir matriz de rigidez global em funcao das barras
def calcular_kg(pontos, elementos):
    """Generate de stiffness matrix.

    Parameters:

    pontos: list of integers
    elementos: list containing elements object type
    """

    assert len(elementos) == len(set(elementos)),\
              "Tem elementos repetidos na lista!"

    kg = np.zeros((ngl*len(pontos), ngl*len(pontos)), dtype=elem_matriz)

    #Acrescenta-se as matrizes de cada elementos a matriz global
    for i in elementos:
        xi, yi, j, k = i.calcular_kge(pontos, kg)

        kg[xi,yi] += i.Kl[j,k]

    return kg

#Funcao para definir a matriz de massas do sistema
def calcular_mg(pontos, elementos):

    assert len(elementos) == len(set(elementos)),\
              "Tem elementos repetidos na lista!"
              
    #zera-se a matriz mg inicialmente
    mg = np.zeros((ngl*len(pontos), ngl*len(pontos)), dtype=elem_matriz)
    #mg = np.zeros((ngl*len(pontos),ngl*len(pontos)), dtype = elem_matriz)

    for i in elementos:
        mg = i.introduzir_me_em_mg(mg)

    return mg
#%% ### MUDADO RECENTEMENTE
def delete_row_lil(matriz0, i):
    # Apaga a linha de uma dada matriz
    mat = matriz0.copy()
    if not isinstance(mat, lil_matrix):
        raise ValueError("works only for LIL format -- use .tolil() first")
    mat.rows = np.delete(mat.rows, i)
    mat.data = np.delete(mat.data, i)
    mat.reshape = (mat.shape[0] - 1, mat.shape[1]) ### MUDADO RECENTEMENTE ATENCAO

    return mat

def remover_linhas_e_colunas(mat, i):
    # Remove as linhas e colunas de indice i de uma dada matriz
    MAT = mat
    MAT = delete_row_lil(MAT, i)
    MAT = MAT.transpose()
    MAT = delete_row_lil(MAT, i)
    MAT = MAT.transpose()

    return MAT
#%%
#Funcao para atribuir as restricoes a matriz de rigidez ou a de forca
def atribuir_restricoes(XX, rest):
    """Remove as linhas e colunas da matriz XX para se implementar as restricoes
    """
    assert len(rest) == len(set(rest))
    rest = list(set(rest))

    if np.shape(XX)[1] > 1 and XX.__class__.__name__ != "lil_matrix":

        XXr = np.delete(XX, rest, 0)
        XXr = np.delete(XXr, rest, 1)

    elif np.shape(XX)[1] == 1:

        XXr = np.delete(XX, rest, 0)

    elif XX.__class__.__name__ == "lil_matrix":
        XXr = XX
        for i in range(len(rest)-1, -1, -1):

            XXr = remover_linhas_e_colunas(XXr, rest[i])
        XXr.tolil()

    else:

        dimensao = len(np.shape(XX))
        print "Dimensao de array invalida!!", dimensao

    return XXr

#Funcao para definir a matriz de forca global
def definir_fg(forcas_nodais, tamanho_pontos): ###############rever nome de funcao
    """Com base numa lista contendo os elementos a seguir esta funcao cria
    o vetor de forcas.

    Convencoes adotadas:
        O numero do ponto deve ser > 0
        Convencao Graus de Liberdade:
            1 = deslocamento na direcao X
            2 = deslocamento na direcao Y
            3 = deslocamento na direcao Z
            4 = rotacao em torno de X
            5 = rotacao em torno de Y
            6 = rotacao em torno de Z

    Parametres:
        forcas_nodais = list
            A lista forcas_nodais deve esta no formato a seguir:
            forcas_nodais = [ [ Numero do ponto de aplicacao da forca,
                                grau de liberdade associado a essa forca,
                                valor desta forca ] ]
        tamanho_pontos : integer
            len( lista_de_pontos )
    """
    fg = np.zeros((ngl*tamanho_pontos, 1), dtype=elem_matriz)

    for i in forcas_nodais:
        fg[(i[0] - 1)* ngl + i[1] - 1, 0] += i[2]

        if i[0] == 0:

            print "Atencao!! Nao existe ponto zero para se aplicar forca!!"

    return fg

#Resolver o sistema de equacoes lineares
def resolver_sistema(A, B):
    # Resolve o sistema de equacoes com uma rotina otimizada da biblioteca numpy

    try:
        resposta = np.linalg.solve(A, B)
    except np.linalg.LinAlgError:
        resposta, residuo, rank, s = np.linalg.lstsq(A, B)

    # Verifica-se se o sistema pode ser mal condicionado
    if np.linalg.cond(A) > 10**15:
        # Para matrizes com elementos de precisao dupla o numero condicional
        # deve ser menor ou igual a 10**15 para garantir respostas adequadas
        # do sistema de equacoes.
        msg = ("ATENCAO! Sistema de equacoes mal condicionado!\n"
               "Condicional da matriz Kgr = %s\n" % np.linalg.cond(A))
        warnings.showwarning(msg, UserWarning, '', -1)
    # Verificacao dos resultados obtidos
    assert np.allclose(np.dot(A, resposta), B)
#    if not np.allclose(np.dot(A,resposta),B):
#        warnings.showwarning(
#           "ATENCAO!! Solucao invalida do sistema de equacoes!\n",\
#           UserWarning,'',-1)

    return resposta
#%%
def encontrar_frequencias(kgr, mgr):
    """
    Calcula a matriz 'v' e o vetor 'b' que soluciona a seguinte equacao
     [ K ] * { v } = { b } * [ m ] * { v }

    K é a matriz de rigidez global com as restrincoes atribuidas
    b é o vetor com o quadrado das frequencias naturais (autovalores)
    m é a matriz global de massas
    v contem os autovetores do sistema
    """
    # Metodo 1 : tratar todas as matrizes como densas
    try:
        mgri = np.linalg.inv(mgr)
    except np.linalg.LinAlgError:
        mgri = np.linalg.pinv(mgr)
    D = np.dot(mgri, kgr)
    b, v = np.linalg.eig(D)

#    # Metodo 2 : converter matris sparse em densa
#    #k = np.shape(kgr)[0]-1
#    mgr = mgr.todense()
#    #minv = np.linalg.pinv(mgr)
#    minv = np.linalg.inv(mgr)
#    A = np.dot(minv,kgr)
#
#    b, v = eigsh(A)
#
#    # Metodo 3 : para resolver matriz sparse
#    k = np.shape(kgr)[0]-1
#
#    b , v = eigsh(kgr,M=mgr,k=k,which="LM")

    w = np.sqrt(b)
    w = np.sort(w)

    return w, v
 #%%
def expandir_xgr(Xr, pontos, rest):
    """Cria o vetor X com os deslocamentos de todos os nos"""

    # Inicia-se X como um vetor de zeros
    X = np.zeros((ngl*len(pontos), 1), dtype=elem_matriz)

    j = 0

    #Acrescenta os termos de Xr em sua respectiva posicao em X
    for i in range(ngl*len(pontos)):

        if i not in rest:

            X[i, 0] = Xr[j, 0]
            j += 1

    return X

def recriar_fg(Kgr, Xr, Fgr, Kg, X):  ### Funcao obsoleta, tavez possa ser removida

    #Calcula a matriz de forcas global do sistema, que agora apresentara as
    #reacoes nos apoios
    Fg = np.dot(Kg, X)

    return Fg
#%%
def calcular_esforcos_barras(barras, X, pontos):

    for elementos in barras:

        i1 = pontos.index(elementos.P1)*ngl #Indice do no 1 na lista de pontos
        i2 = pontos.index(elementos.P2)*ngl #Indice do no 2 na lista de pontos

        elementos.Xge = np.zeros(2*ngl, dtype=elem_matriz)

        for i in range(2*ngl):

            if i < ngl:

                elementos.Xge[i] = X[i1 + i]

            else:

                elementos.Xge[i] = X[i2 + (i-ngl)]

        Fl = np.dot(elementos.Ke, elementos.T)
        Fl = np.dot(Fl, elementos.Xge)

        elementos.Fl = Fl

def calcular_esforcos_CST(lista_elementos, X, pontos):

    for elementos in lista_elementos:

        elementos.Xge = np.zeros(elementos.npontos*ngl, dtype=elem_matriz)

        for i in range(elementos.npontos*ngl):

            j1 = (i - i % ngl)/ngl
            x = elementos.ref_pontos[j1]+i%ngl

            elementos.Xge[i] = X[x]

        Fl = np.dot(elementos.D, elementos.B)
        Fl = np.dot(Fl, elementos.Xge)

        elementos.Fl = Fl
#%%
def metodo_da_diferenca_central(m, k, f, tf, dt):

    assert isinstance(tf, float) and isinstance(dt, float), \
                     "tf e dt devem ser numeros reais!"

    # Calculo do numero de passos a serem executados
    n_dt = int(tf/dt)

    # Cria-se a lista que recebera os vetores com os deslocamentos
    x = []
    x.append(np.zeros((np.shape(m)[0], 1)), dtype=elem_matriz)
    v0 = np.zeros((np.shape(m)[0], 1), dtype=elem_matriz)

    # Calculo da matriz de massa inversa
    m_inv = np.linalg.inv(m)
    m_inv_k = np.dot(m_inv, k)

    # Primeira iteracao
    aux1 = ((dt**2.0)/2.0)* np.dot(m_inv, f[0])
    aux2 = (((dt**2.0)/2.0)* np.dot(m_inv_k, x[0]))
    aux3 = (dt* v0 + x[0])
    aux = aux1 - aux2 + aux3

    x.append(aux)

    t = [0.0, dt]

    for i in range(1, n_dt):
        aux1 = dt**2.0 * np.dot(m_inv, f[i])
        aux2 = dt**2.0 * np.dot(m_inv_k, x[i])
        aux3 = 2 * x[i] - x[i-1]
        aux = aux1 - aux2 + aux3
        x.append(aux)
        t.append((i+1)*dt)

    return x, t

def metodo_de_newmark(m, k, f, tf, dt, beta=0.25, gama=0.5):
    """
    Seguindo o procedimento apresentado por Logan(2007) capitulo 16:
    Structural Dynamics and Time-Dependent Heat Transfer

    Diferenca central:
        beta = 0  #Erro pois beta tem que ser diferente de zero
        gama = 0.5

    Metodo da aceleracao linear:
        beta = 1/6
        gama = 0.5

    Metodo da aceleracao media:
        beta = 0.25
        gama = 0.5
    """
    c = np.zeros((np.shape(m)[0], np.shape(m)[0]),
                 dtype=elem_matriz)
    #Inserir esse termo na formulacao depois


    assert isinstance(tf, float) and isinstance(dt, float), \
                     "tf e dt devem ser numeros reais!"

    # Calculo da matriz de massa inversa
    try:
        m_inv = np.linalg.inv(m)
    except np.linalg.LinAlgError:
        m_inv = np.linalg.pinv(m)

    # Calculo do numero de passos a serem executados
    n_dt = int(tf/dt)

    # 1o PASSO
    # Cria-se a lista que recebera os vetores com os deslocamentos
    x = []
    x.append(np.zeros((np.shape(m)[0], 1), dtype=elem_matriz))
    v = []
    v.append(np.zeros((np.shape(m)[0], 1), dtype=elem_matriz))
    a = []
    t = [0.0]

    # 2o PASSO
    a.append(np.dot(m_inv, f[0] - np.dot(c, v[0]) - np.dot(k, x[0])))

    aux = k + (1.0/(beta*dt**2.0))*m + (gama/(beta*dt)) * c
    try:
        k1_inv = np.linalg.inv(aux)
    except np.linalg.LinAlgError:
        k1_inv = np.linalg.pinv(aux)

    for i in range(n_dt):

        # 3o PASSO

        aux = x[i] + dt*v[i] + (0.5 - beta)*(dt**2.0)*a[i]
        f1 = f[i+1] + (1.0/(beta*dt**2.0)) * np.dot(m, aux)

        aux = np.dot(k1_inv, f1)

        x.append(aux)

        # 4o PASSO

        aux = x[i+1] - x[i] - dt*v[i] - dt**2.0 * (0.5 - beta)*a[i]
        aux = (1.0/(beta*dt**2.0)) * aux

        a.append(aux)

        # 5o PASSO

        aux = v[i] + dt * ((1.0 - gama) * a[i] + gama * a[i+1])
        v.append(aux)

        t.append(dt*(i+1))

    return x, t

def f_dinamico_func(f, w, t, dt, func):
    """
    Funcao para criar lista onde cada indice corresponde a o vetor forca em um
    instante de tempo.
    Note que 'func' eh uma funcao que pode ser 'sin', 'cos', etc..
    """

    f_dinamico = []

    for i in range(int(t/dt)+1):

        f_dinamico.append(f * func(w*i*dt))

    return f_dinamico

def resposta_din_ponto_i(num_linha, x_din):

    resposta = []

    for i in range(len(x_din)):

        resposta.append(x_din[i][num_linha])

    return resposta
#%%
def gerar_matrizes(pontos, barras, restricoes, forcas, massas_pontuais=[]):
    """Gera as matrizes de rigidez, massa e de forca do sistema.

    Args:
        pontos (list) : lista com as coordenadas dos pontos
        barras (list) : lista com os objetos do tipo Barra
        restricoes (list) : lista com os valores de linhas e colunas a serem \
                            apagadas das matrizes geradas
        forcas (list) : lista referente os esforcos solicitantes nodais que\
                        devem está como a seguir.
            forcas_nodais = [ [ Numero do ponto de aplicacao da forca,
                                grau de liberdade associado a essa forca,
                                valor desta forca ] ]
        massas_pontuais (lista) : lista com as massas pontuais.

    Returns:
        kg,mg,fg (np.arrays) : Matrizes de rigidez, massa e vetor de forças
        kgr,mgr,fgr (np.arrays) : Matrizes de rigidez, massa e vetor de forças
                                  com os graus de liberdades restritos removidos
    """
    kg = calcular_kg(pontos, barras)

    kgr = atribuir_restricoes(kg, restricoes)

    mg = calcular_mg(pontos, barras)

    for i in massas_pontuais:

        mg[i[1], i[1]] = mg[i[1], i[1]] + i[0]

    mgr = atribuir_restricoes(mg, restricoes)

    fg = definir_fg(forcas, len(pontos))

    fgr = atribuir_restricoes(fg, restricoes)

    return kg, kgr, mg, mgr, fgr, fg

def analise_estatica(kg, kgr, fgr, pontos, restricoes, forcas):

    xgr = resolver_sistema(kgr, fgr)

    x = expandir_xgr(xgr, pontos, restricoes)

    fg_c_reacoes = recriar_fg(kgr, xgr, fgr, kg, x)

    return x, fg_c_reacoes

def analise_dinamica(kgr, mgr, fgr, freq_carreg, t_carreg, dt,
                     funcao, pontos, restricoes):

    w, v = encontrar_frequencias(kgr, mgr)

    forcas_dinamicas = f_dinamico_func(
        fgr, freq_carreg, t_carreg, dt, funcao)

    x_din, t = metodo_de_newmark(mgr, kgr, forcas_dinamicas, t_carreg, dt)

    x_din2 = []
    for i in range(len(x_din)):

        x_din2.append(expandir_xgr(x_din[i], pontos, restricoes))

    return w, v, x_din2, t, forcas_dinamicas
#%%
if __name__ == "__main__":

    # Exemplo Lista 05 - dinamica das estruturas
    from math import cos
    k = 12.0*28*10**9*0.4*0.2**3/12.0/3**3
    K = np.matrix(([8.0*k, -4.0*k], [-4.0*k, 4.0*k]))
    M = np.eye(2)*6750
    f_din = f_dinamico_func(np.matrix(([[1000.0], [1000.0]])),
                            10.0, 1.0, 0.001, cos)
    #x_din,t = metodo_da_diferenca_central(M,K,f_din,1.0,0.001)
    x_din, t = metodo_de_newmark(M, K, f_din, 1.0, 0.001)
    testando = resposta_din_ponto_i(1, x_din)
    from matplotlib.pyplot import plot
    plot(t, testando)
