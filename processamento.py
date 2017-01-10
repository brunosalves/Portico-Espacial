# -*- coding: utf-8 -*-
"""
UNIVERSIDADE FEDERAL DE PERNAMBUCO
Autor: Bruno Sampaio Alves
Disciplina: Metodo dos elementos Finitos
Professor: Paulo Marcelo Vieira Ribeiro

                        FUNCOES DE PROCESSAMENTO
"""
from settings import ngl,ngls,elem_matriz
from scipy.sparse import lil_matrix
#from scipy.sparse.linalg import eigsh, inv
#from scipy.linalg import solve_banded # Para solucionar um problema com sparse_matrix
import numpy as np

#Funcao para definir matriz de rigidez global em funcao das barras
def calcular_kg(pontos, elementos):    
    """Generate de stiffness matrix.
    
    Parameters:
    
    pontos: list of integers
    elementos: list containing elements object type
    
    """
    
    Kg = np.zeros((ngl*len(pontos),ngl*len(pontos)), dtype = elem_matriz)    
    
    #Acrescenta-se as matrizes de cada elementos a matriz global
    for i in elementos:
        Kg = i.calcular_kge(pontos,Kg)
        
    return Kg

#Funcao para definir a matriz de massas do sistema
def calcular_mg(pontos, elementos):
    
    #zera-se a matriz mg inicialmente
    mg = np.zeros((ngl*len(pontos),ngl*len(pontos)), dtype = elem_matriz)    
    #mg = np.zeros((ngl*len(pontos),ngl*len(pontos)), dtype = elem_matriz)    
        
    
    for i in elementos:
 
        mg = i.introduzir_me_em_mg(mg)

    return mg

def delete_row_lil(m, i):
    mat = m
    if not isinstance(mat, lil_matrix):
        raise ValueError("works only for LIL format -- use .tolil() first")
    mat.rows = np.delete(mat.rows, i)
    mat.data = np.delete(mat.data, i)
    mat._shape = (mat._shape[0] - 1, mat._shape[1])
    
    return mat

def remover_linhas_e_colunas(mat,i):
      
    m = mat
    m = delete_row_lil(m, i)
    m = m.transpose()
    m = delete_row_lil(m, i)
    m = m.transpose()
    
    return m
  
#Funcao para atribuir as restricoes a matriz de rigidez ou a de forca    
def atribuir_restricoes(XX,rest):

    rest = list(set(rest))
    """
    Remove as linhas e colunas da matriz XX para se implementar as restricoes
    """
    if np.shape(XX)[1] > 1 and XX.__class__.__name__ != "lil_matrix":
        
        XXr = np.delete(XX, rest, 0)
        XXr = np.delete(XXr, rest, 1)
        
    elif np.shape(XX)[1] == 1:
        
        XXr = np.delete(XX, rest, 0)
    
    elif XX.__class__.__name__ == "lil_matrix":    
        XXr = XX
        for i in range(len(rest)-1,-1,-1):
            
            XXr = remover_linhas_e_colunas(XXr,rest[i])
        XXr.tolil()

    else:
        
        dimensao = len(np.shape(XX))        
        print "Dimensao de array invalida!!", dimensao
    
    return XXr

#Funcao para definir a matriz de forca global
def definir_fg(forcas_nodais,tamanho_pontos): ###############rever nome de funcao

    fg = np.zeros((ngl*tamanho_pontos,1), dtype = elem_matriz)
    
    for i in forcas_nodais:
        fg[ ( i[0] - 1 ) * ngl + i[1] - 1 , 0 ] += i[2]
        
        if i[0] == 0:

            print "Atencao!! Nao existe ponto zero para se aplicar forca!!"
       
    return fg
    
#Resolver o sistema de equacoes lineares
def resolver_sistema(A,B):
    
    # Resolve o sistema de equacoes com uma rotina otimizada da biblioteca numpy   
    
    try:
        resposta = np.linalg.solve(A, B)
    except np.linalg.LinAlgError:
        A2 = np.linalg.pinv(A)
        resposta = np.dot(A2,B)        
    
    return resposta

def encontrar_frequencias(kgr,mgr):
    
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
    D = np.dot(mgri ,kgr)
    b,v = np.linalg.eig(D)    
    
    """
    # Metodo 2 : converter matris sparse em densa
    #k = np.shape(kgr)[0]-1
    mgr = mgr.todense()
    #minv = np.linalg.pinv(mgr)
    minv = np.linalg.inv(mgr)
    A = np.dot(minv,kgr)
    
    b, v = eigsh(A)   
    """
    """
    # Metodo 3 : para resolver matriz sparse
    k = np.shape(kgr)[0]-1
    
    b , v = eigsh(kgr,M=mgr,k=k,which="LM")
    """
    w = np.sqrt(b)
    w = np.sort(w)
    return w, v
    
def expandir_xgr(Xr,pontos, rest):
    """Cria o vetor X com os deslocamentos de todos os nos"""
    
    # Inicia-se X como um vetor de zeros
    X = np.zeros((ngl*len(pontos),1), dtype = elem_matriz)
    
    j = 0
    
    #Acrescenta os termos de Xr em sua respectiva posicao em X
    for i in range(ngl*len(pontos)):
        
        if i not in rest:
            
            X[i,0] = Xr[j,0]
            j = j + 1    
            
    return X

def recriar_fg(Kgr,Xr,Fgr,Kg,X):  ### Verificar nome de funcao muito similar
    
    # Verifica se os valores de Xr encontrados multiplicados por Kgr dao Fgr
    if np.allclose(np.dot(Kgr,Xr),Fgr) == False:
        print "ERRO!!! Solucao invalida do sistema de equacoes!"
    
    #Calcula a matriz de forcas global do sistema, que agora apresentara as 
    #reacoes nos apoios    
    Fg =  np.dot(Kg,X)  
    
    return Fg    

def calcular_esforcos_barras(barras,X,pontos):   
    
    for elementos in barras:
        
        i1 = pontos.index(elementos.P1)*ngl #Indice do no 1 na lista de pontos
        i2 = pontos.index(elementos.P2)*ngl #Indice do no 2 na lista de pontos        
        
        elementos.Xge = np.zeros(2*ngl, dtype = elem_matriz)
        
        for i in range(2*ngl):
            
            if i < ngl:
                
                elementos.Xge[i] = X[i1 + i]
                
            else:
                
                elementos.Xge[i] = X[i2 + (i-ngl)]
                
        Fl = np.dot(elementos.Ke,elementos.T)
        Fl = np.dot(Fl,elementos.Xge)
             
        elementos.Fl = Fl
             
def calcular_esforcos_CST(lista_elementos,X,pontos):   
    
    for elementos in lista_elementos:
        
        elementos.Xge = np.zeros(elementos.npontos*ngl, dtype = elem_matriz)
        
        for i in range(elementos.npontos*ngl):
            
            j1 = (i-i%ngl)/ngl
            x = elementos.ref_pontos[j1]+i%ngl            
            
            elementos.Xge[i] = X[x]
              
        Fl = np.dot(elementos.D,elementos.B)
        Fl = np.dot(Fl,elementos.Xge)
            
        elementos.Fl = Fl
        
def metodo_da_diferenca_central(m,k,f,tf,dt):
    
    tf = float(tf)
    dt = float(dt)    
    
    # Calculo do numero de passos a serem executados
    n_dt = int(tf/dt)
    
    # Cria-se a lista que recebera os vetores com os deslocamentos
    x = []
    x.append(np.zeros((np.shape(m)[0],1)), dtype = elem_matriz)
    v0 = np.zeros((np.shape(m)[0],1), dtype = elem_matriz)
    
    # Calculo da matriz de massa inversa
    m_inv = np.linalg.inv(m)
    m_inv_k = np.dot(m_inv,k)
    
    # Primeira iteracao
    aux1 = ((dt**2.0)/2.0) * np.dot(m_inv,f[0])
    aux2 = (((dt**2.0)/2.0) * np.dot( m_inv_k, x[0] ))
    aux3 = ( dt* v0 + x[0] )
    aux = aux1 - aux2 + aux3
         
    x.append(aux)

    t = [0.0,dt]    
    
    for i in range(1,n_dt):
        aux1 = dt**2.0 * np.dot(m_inv,f[i])  
        aux2 = dt**2.0 * np.dot( m_inv_k, x[i] )
        aux3 = 2 * x[i] - x[i-1]
        aux = aux1 - aux2 + aux3
        x.append(aux)
        t.append((i+1)*dt)
    
    return x,t

def metodo_de_newmark(m,k,f,tf,dt,beta=0.25,gama=0.5):
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
    c=np.zeros((np.shape(m)[0],np.shape(m)[0]), dtype = elem_matriz) #Inserir esse termo na formulacao depois
    tf = float(tf)
    dt = float(dt)   
    
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
    x.append(np.zeros((np.shape(m)[0],1), dtype = elem_matriz))
    v = []    
    v.append(np.zeros((np.shape(m)[0],1), dtype = elem_matriz))
    a = []
    t = [0.0]
    
    # 2o PASSO    
    a.append( np.dot( m_inv, f[0] - np.dot(c,v[0]) - np.dot(k,x[0])) )
    
    aux = k + (1.0/(beta*dt**2.0))*m + (gama/(beta*dt)) * c
    try:    
        k1_inv = np.linalg.inv( aux )   
    except np.linalg.LinAlgError:  
        k1_inv = np.linalg.pinv( aux )   
        
    for i in range(n_dt):
        
        # 3o PASSO
               
        aux = x[i] + dt*v[i] + (0.5 - beta)*(dt**2.0)*a[i]
        f1 = f[i+1] + (1.0/(beta*dt**2.0)) * np.dot(m, aux)
        
        aux = np.dot( k1_inv , f1 )
        
        x.append(aux)
        
        # 4o PASSO
    
        aux = x[i+1] - x[i] - dt*v[i] - dt**2.0 * (0.5 - beta)*a[i]
        aux = (1.0/(beta*dt**2.0)) * aux
        
        a.append( aux )
        
        # 5o PASSO
        
        aux = v[i] + dt * ( (1.0 - gama) * a[i] + gama * a[i+1] )
        v.append(aux)
        
        t.append(dt*(i+1))
    
    return x,t
    
def f_dinamico_func(f,w,t,dt,func):
    """
    Funcao para criar lista onde cada indice corresponde a o vetor forca em um
    instante de tempo.
    Note que 'func' eh uma funcao que pode ser 'sin', 'cos', etc..
    """
        
    f_dinamico = []
    
    for i in range(int(t/dt)+1):
        
        f_dinamico.append(f * func(w*i*dt))
        
    return f_dinamico
    
def resposta_din_ponto_i(num_linha,x_din):
    
    resposta = []
    
    for i in range(len(x_din)):
        
        resposta.append(x_din[i][num_linha])
        
    return resposta

if __name__ == "__main__":
    
    # Exemplo Lista 05 - dinamica das estruturas
    from math import cos
    k=12.0*28*10**9*0.4*0.2**3/12.0/3**3
    K=np.matrix(([8.0*k,-4.0*k],[-4.0*k,4.0*k]))
    M=np.eye(2)*6750
    f_din = f_dinamico_func(np.matrix(([[1000.0],[1000.0]])),10.0,1.0,0.001,cos)
    #x_din,t = metodo_da_diferenca_central(M,K,f_din,1.0,0.001)
    x_din,t = metodo_de_newmark(M,K,f_din,1.0,0.001)
    testando = resposta_din_ponto_i(1,x_din)    
    from matplotlib.pyplot import plot
    plot(t,testando)