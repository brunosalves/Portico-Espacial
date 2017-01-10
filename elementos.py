# -*- coding: utf-8 -*-
"""
UNIVERSIDADE FEDERAL DE PERNAMBUCO
Autor: Bruno Sampaio Alves
Disciplina: Metodo dos Elementos Finitos
Professor: Paulo Marcelo Vieira Ribeiro

            Tipos de elementos
"""
from settings import ngl, elem_matriz
from sympy import symbols, diff, Matrix, integrate
from scipy.sparse import lil_matrix

import numpy as np

##### Elementos de barras tridimensionais

class Barra(object):

    npontos = 2

    def __init__(self, P1, P2, A, E, Iz, Iy, G, J, p):
        self.P1 = P1
        self.P2 = P2
        self.lista_pontos = [P1,P2]
        self.A = A
        self.E = E
        self.Iz = Iz
        self.Iy = Iy
        self.G = G
        self.J = J
        self.L = np.sqrt(((np.array(P2)-np.array(P1))*\
        (np.array(P2)-np.array(P1))).sum()) #Modulo do vetor dado por P2-P1
        self.p = p #densidade por metro linear da estrutura        
        
        #calcula-se a inclinacao em relacao a x
        self.l = (P2[0] - P1[0]) / self.L
        self.m = (P2[1] - P1[1]) / self.L
        self.n = (P2[2] - P1[2]) / self.L
        
        self.Ke = self.calcular_ke()
        self.Lambda = self.calcular_lambda()        
        self.T = self.calcular_T()
        self.me = self.calcular_me()

    def definir_ref_pontos(self,pontos):
        
        ref_pontos = []
        
        for i in self.lista_pontos:            
            ref_pontos.append(pontos.index(i)*ngl)
            
        return ref_pontos

    def calcular_ke(self):
        #Constante que sera multiplicada pela matriz de rigidez adimensional
        EA_L = self.E * self.A / self.L
        EIz_L3 = self.E * self.Iz / ( self.L ** 3 )       
        EIy_L3 = self.E * self.Iy / ( self.L ** 3)        
        GJ_L = self.G * self.J / self.L
        
        
        L = self.L
        
        #Matriz de rigidez
        Ke = np.zeros( (12,12) , dtype = elem_matriz)       
        
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

    def calcular_lambda(self):
        
        D = ( self.l ** 2 + self.m ** 2 ) ** 0.5        

        Lambda = np.zeros( ( 3 , 3 ))
        Lambda[0,0] = self.l
        Lambda[0,1] = self.m
        Lambda[0,2] = self.n

        if self.l == 0 and self.m == 0:
            Lambda[1,1] = 1
            Lambda[1,0] = 0
            Lambda[2,0] = (-1.0)*self.n
            Lambda[2,1] = 0
        else:
            Lambda[1,1] = self.l/D
            Lambda[1,0] = -self.m/D
            Lambda[2,0] = -self.l*self.n/D
            Lambda[2,1] = -self.m*self.n/D
        
        Lambda[1,2] = 0
        Lambda[2,2] = D      
        
        return Lambda

    def calcular_T(self):
        
        #Matriz rotacao
        T = np.zeros((12,12))

        Lambda = self.Lambda        
        
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
        
    def calcular_kge(self,pontos,Kge):         
        """Calcula a matriz kl do elemento e introduz na matriz Kge do sistema
        """        
        
           
        # [Kl] =  [T]^t . [Ke] . [T]   
        Tt = np.matrix.transpose(self.T)
        kl = np.dot(Tt,self.Ke)
        kl = np.dot(kl,self.T) #Matriz local rotacionada para os eixos globais    
        self.Kl = kl    
        
        #definir a matriz global so com a parte deste elemento
        
        ref_pontos = self.definir_ref_pontos(pontos)
        self.ref_pontos = ref_pontos
        
        #Introduzir a KL em Kge
        for i in range(len(ref_pontos)*ngl):
               
            for j in range(len(ref_pontos)*ngl):                   
                
                j1 = (i-i%ngl)/ngl
                k1 = (j-j%ngl)/ngl
                
                x = ref_pontos[j1]+i%ngl
                y = ref_pontos[k1]+j%ngl
                
                Kge[x,y] =+ kl[i,j]
                        
        return Kge

    def calcular_me(self):
        
        me = lil_matrix((12,12))

        if self.A != 0:
            # Diagonal Principal        
            me[0,0] = 140.0
            me[1,1] = 156.0
            me[2,2] = me[1,1]
            me[3,3] = 140.0*self.J/self.A
            me[4,4] = 4.0 * self.L**2.0
            me[5,5] = me[4,4]
            me[6,6] = me[0,0]
            me[7,7] = me[1,1]
            me[8,8] = me[1,1]
            me[9,9] = me[3,3]
            me[10,10] = me[4,4]
            me[11,11] = me[4,4]
            
            # Triangulo superior
            me[0,6] = 70.0
            me[1,5] = 22.0*self.L
            me[1,7] = 54.0
            me[1,11] = -13.0*self.L
            me[2,4] = (-1.0) * me[1,5]
            me[2,8] = me[1,7]
            me[2,10] = (-1) * me[1,11]
            me[3,9] = 70.0 * self.J / self.A
            me[4,8] = me[1,11]
            me[4,10] = -3.0 * self.L**2.0
            me[5,7] = (-1.0) * me[1,11]
            me[5,11] = me[4,10]
            me[7,11] = (-1.0) * me[1,5]
            me[8,10] = me[1,5]
            
            me[6,0] = me[0,6]
            me[5,1] = me[1,5]
            me[7,1] = me[1,7]
            me[11,1] = me[1,11]
            me[4,2] = me[2,4]
            me[8,2] = me[2,8]
            me[10,2] = me[2,10]
            me[9,3] = me[3,9]
            me[8,4] = me[4,8]
            me[10,4] = me[4,10]
            me[7,5] = me[5,7]
            me[11,5] = me[5,11]
            me[11,7] = me[7,11]
            me[10,8] = me[8,10]
    
            me = me * self.p * self.L / 420.0
         
        return me

    def introduzir_me_em_mg(self,mg):
        
        #Introduzir a me em mge

        for i in range(len(self.ref_pontos)*ngl):
               
            for j in range(len(self.ref_pontos)*ngl):                   
                
                j1 = (i-i%ngl)/ngl
                k1 = (j-j%ngl)/ngl
                
                x = self.ref_pontos[j1]+i%ngl
                y = self.ref_pontos[k1]+j%ngl
                
                mg[x,y] = mg[x,y] + self.me[i,j]                    

        return mg

############## Definir elementos CST bidimensionais ##################
class ElementosPlanos(object):

    def calcular_D(self):
        #Para o Estado plano de tensoes        
        D = np.zeros((3,3))
        D[0,0] = 1
        D[0,1] = self.ni
        D[1,0] = self.ni
        D[1,1] = 1
        D[2,2] = (1 - self.ni) / 2
        D = self.E/(1-self.ni**2)*D
           
        """ #Para o Estado plano de deformacoes
        D = np.zeros((3,3))
        D[0,0] = 1 - self.ni
        D[0,1] = self.ni
        D[1,0] = self.ni
        D[1,1] = 1 - self.ni
        D[2,2] = (1 - 2 * self.ni) / 2
        D = self.E / ((1+self.ni)*(1-2*self.ni))
        """
        
        return D

    def definir_ref_pontos(self,pontos):
        
        ref_pontos = []
        
        for i in self.lista_pontos:            
            ref_pontos.append(pontos.index(i)*ngl)
            
        return ref_pontos

class ElementoCST(ElementosPlanos):

    npontos = 3

    def __init__(self, P1, P2, P3, E, ni, t): 
        
        self.P1 = P1
        self.P2 = P2
        self.P3 = P3
        self.lista_pontos = [P1,P2,P3]
        self.E = E
        self.ni = ni
        self.t = t
        self.A = self.calcular_A()
        self.B = self.calcular_B()
        self.D = self.calcular_D()
        
    def calcular_A(self):
        #Define-se A como uma matriz de zeros
        A = np.zeros((3,3))
        
        #Primeira coluna
        A[0,0] = 1
        A[1,0] = 1
        A[2,0] = 1
        
        #Segunda coluna
        A[0,1] = self.P1[0]
        A[1,1] = self.P2[0]
        A[2,1] = self.P3[0]
        
        #Terceira coluna
        A[0,2] = self.P1[1]
        A[1,2] = self.P2[1]
        A[2,2] = self.P3[1]
        
        A = 0.5*np.linalg.det(A)
        #Atribui-se A ao elemento CST
        
        return A
        
    def calcular_B(self):
        #Define-se B como uma matriz de zeros
        B = np.zeros((3,6))

        #Primeira linha de B
        B[0,0] = self.P2[1]-self.P3[1]
        B[0,2] = self.P3[1]-self.P1[1]
        B[0,4] = self.P1[1]-self.P2[1]
        
        #Segunda linha de B
        B[1,1] = self.P3[0]-self.P2[0]
        B[1,3] = self.P1[0]-self.P3[0]
        B[1,5] = self.P2[0]-self.P1[0]
        
        #Terceira linha
        B[2,0] = self.P3[0]-self.P2[0]
        B[2,1] = self.P2[1]-self.P3[1]
        B[2,2] = self.P1[0]-self.P3[0]
        B[2,3] = self.P3[1]-self.P1[1]
        B[2,4] = self.P2[0]-self.P1[0]
        B[2,5] = self.P1[1]-self.P2[1]
        
        #Atribui-se o B ao elemento CST
        B = ( 1.0 / ( 2 * self.A ) ) * B
        
        return B

    def calcular_kge(self, pontos, Kge):
        
        constante = self.t * self.A
        
        kl = np.dot(np.transpose(self.B),self.D)
        kl = np.dot(kl,self.B)
        kl = constante*kl
        
        self.Kl = kl
        
        ref_pontos = self.definir_ref_pontos(pontos)
        self.ref_pontos = ref_pontos
        
        #Introduzir a KL em Kge
        for i in range(len(ref_pontos)*ngl):
               
            for j in range(len(ref_pontos)*ngl):                   
                
                j1 = (i-i%ngl)/ngl
                k1 = (j-j%ngl)/ngl
                
                x = ref_pontos[j1]+i%ngl
                y = ref_pontos[k1]+j%ngl
                
                Kge[x,y] = Kge[x,y] + kl[i,j]                    

        return Kge
        
class ElementoQ4(ElementosPlanos):

    npontos = 4

    def __init__(self, P1, P2, P3, P4, E, ni, t): 
        
        self.P1 = P1
        self.P2 = P2
        self.P3 = P3
        self.P4 = P4
        self.lista_pontos = [P1,P2,P3,P4]
        self.E = E
        self.ni = ni
        self.t = t
        self.D = self.calcular_D()
        self.kl = self.calcular_kl()
    
    
    def calcular_kl(self):
        
        t =  self.t
        x1,y1 = self.P1[0],self.P1[1]
        x2,y2 = self.P2[0],self.P2[1]
        x3,y3 = self.P3[0],self.P3[1]
        x4,y4 = self.P4[0],self.P4[1]

        xmin = min(x1,x2,x3,x4)
        xmax = max(x1,x2,x3,x4)
        ymin = min(y1,y2,y3,y4)
        ymax = max(y1,y2,y3,y4)

        u1,u2,u3,u4 = symbols('u1 u2 u3 u4')
        x, y = symbols('x y')
        
        M= Matrix([ [1, x1, y1, x1*y1],
                    [1, x2, y2, x2*y2],
                    [1, x3, y3, x3*y3],
                    [1, x4, y4, x4*y4] ])
        
        C = M**-1 * Matrix([[u1],[u2],[u3],[u4]])
        
        U = Matrix([1, x, y, x*y])        
        
        u = U.T * C
        u = u[0]        
        
        N1 = diff(u,u1)
        N2 = diff(u,u2)
        N3 = diff(u,u3)
        N4 = diff(u,u4)
        
        B = Matrix([ [diff(N1,x),          0, diff(N2,x),          0, diff(N3,x),          0, diff(N4,x),          0],
                     [         0, diff(N1,y),          0, diff(N2,y),          0, diff(N3,y),          0, diff(N4,y)],
                     [diff(N1,y), diff(N1,x), diff(N2,y), diff(N2,x), diff(N3,y), diff(N3,x), diff(N4,y), diff(N4,x)]
        
        ])

        D = self.D                
        
        return t * integrate( B.T * D * B, (x,xmin,xmax) , (y,ymin,ymax) )

    
    def calcular_kge(self,pontos,Kge):
        
        kl = self.kl
        
        ref_pontos = self.definir_ref_pontos(pontos)
        self.ref_pontos = ref_pontos
        
        #Introduzir a KL em Kge
        for i in range(len(ref_pontos)*ngl):
               
            for j in range(len(ref_pontos)*ngl):                   
                
                j1 = (i-i%ngl)/ngl
                k1 = (j-j%ngl)/ngl
                
                x = ref_pontos[j1]+i%ngl
                y = ref_pontos[k1]+j%ngl
                
                Kge[x,y] = Kge[x,y] + kl[i,j]                    

        return Kge
        

class ElementoLST(ElementosPlanos):
    
    npontos = 6

    def __init__(self, P1, P2, P3, P4, P5, P6, E, ni, t): 
        
        self.P1 = P1
        self.P2 = P2
        self.P3 = P3
        self.P4 = P4
        self.P5 = P5
        self.P6 = P6
        self.lista_pontos = [P1,P2,P3,P4,P5,P6]
        self.E = E
        self.ni = ni
        self.t = t
        self.D = self.calcular_D()
        self.kl = self.calcular_kl()
    
    
    def calcular_kl(self):
        
        t =  self.t
        x1,y1 = self.P1[0],self.P1[1]
        x2,y2 = self.P2[0],self.P2[1]
        x3,y3 = self.P3[0],self.P3[1]
        x4,y4 = self.P4[0],self.P4[1]
        x5,y5 = self.P5[0],self.P5[1]
        x6,y6 = self.P6[0],self.P6[1]

        u1,u2,u3,u4,u5,u6 = symbols('u1 u2 u3 u4 u5 u6')
        x, y = symbols('x y')
        
        M= Matrix([ [1, x1, y1, x1**2.0, x1*y1, y1**2.0],
                    [1, x2, y2, x2**2.0, x2*y2, y2**2.0],
                    [1, x3, y3, x3**2.0, x3*y3, y3**2.0],
                    [1, x4, y4, x4**2.0, x4*y4, y4**2.0],
                    [1, x5, y5, x5**2.0, x5*y5, y5**2.0],
                    [1, x6, y6, x6**2.0, x6*y6, y6**2.0] ])
        
        C = M**-1 * Matrix([[u1],[u2],[u3],[u4],[u5],[u6]])
        
        U = Matrix([1, x, y, x**2.0, x*y, y**2.0])        
        
        u = U.T * C
        u = u[0]
        
        N1 = diff(u,u1)
        N2 = diff(u,u2)
        N3 = diff(u,u3)
        N4 = diff(u,u4)
        N5 = diff(u,u5)
        N6 = diff(u,u6)
        
        B = Matrix([ [diff(N1,x),          0, diff(N2,x),          0, diff(N3,x),          0, diff(N4,x),          0, diff(N5,x),          0, diff(N6,x),          0],
                     [         0, diff(N1,y),          0, diff(N2,y),          0, diff(N3,y),          0, diff(N4,y),          0, diff(N5,y),          0, diff(N6,y)],
                     [diff(N1,y), diff(N1,x), diff(N2,y), diff(N2,x), diff(N3,y), diff(N3,x), diff(N4,y), diff(N4,x), diff(N5,y), diff(N5,x), diff(N6,y), diff(N6,x)]
        
        ])
        
        D = self.D

        if y2 > y1:

            #a = y1 + (x - x1) * (y2 - y1)/(x2 - x1)
            #b = y1 + (x - x1) * (y3 - y1)/(x3 - x1)
            c = y2 + (x - x2) * (y3 - y2)/(x3 - x2)

            integral =  t * integrate( B.T * D * B, (y,y1,c), (x,x3,x1) )

        if y2 == y1:

            #a = y1 + (x - x1) * (y2 - y1)/(x2 - x1)
            b = y1 + (x - x1) * (y3 - y1)/(x3 - x1)
            #c = y2 + (x - x2) * (y3 - y2)/(x3 - x2)
            
            integral = t * integrate( B.T * D * B, (y,b,y1), (x,x3,x1) )        
        
        return integral

    
    def calcular_kge(self,pontos,Kge):
        
        kl = self.kl
        
        ref_pontos = self.definir_ref_pontos(pontos)
        self.ref_pontos = ref_pontos
        
        #Introduzir a KL em Kge
        for i in range(len(ref_pontos)*ngl):
               
            for j in range(len(ref_pontos)*ngl):                   
                
                j1 = (i-i%ngl)/ngl
                k1 = (j-j%ngl)/ngl
                
                x = ref_pontos[j1]+i%ngl
                y = ref_pontos[k1]+j%ngl
                
                Kge[x,y] = Kge[x,y] + kl[i,j]                    

        return Kge
        
        
class ElementoQ8(ElementosPlanos):

    npontos = 8

    def __init__(self, P1, P2, P3, P4, P5, P6, P7, P8, E, ni, t): 
        
        self.P1 = P1
        self.P2 = P2
        self.P3 = P3
        self.P4 = P4
        self.P5 = P5
        self.P6 = P6
        self.P7 = P7
        self.P8 = P8
        self.lista_pontos = [P1,P2,P3,P4,P5,P6,P7,P8]
        self.E = E
        self.ni = ni
        self.t = t
        self.D = self.calcular_D()
        self.kl = self.calcular_kl()
    
    
    def calcular_kl(self):
        
        t =  self.t

        x5,x7 = self.P5[0],self.P7[0]
        y6,y8 = self.P6[1],self.P8[1]    

        epslon, ni = symbols('epslon ni')
        
        N1 = (1.0/4.0)*(1.0-epslon)*(ni-1.0)*(epslon+ni+1.0)
        N2 = (1.0/4.0)*(1.0+epslon)*(ni-1.0)*(-epslon+ni+1.0)
        N3 = (1.0/4.0)*(1.0+epslon)*(ni+1.0)*(epslon+ni-1.0)
        N4 = (1.0/4.0)*(-1.0+epslon)*(ni+1.0)*(epslon-ni+1.0)
        N5 = (0.5)*(1.0-ni)*(1.0-epslon**2.0)
        N6 = (0.5)*(1.0+epslon)*(1.0-ni**2.0)
        N7 = (0.5)*(1.0+ni)*(1.0-epslon**2.0)
        N8 = (0.5)*(1.0-epslon)*(1.0-ni**2.0)
        
        B = Matrix([ [diff(N1,epslon),               0, diff(N2,epslon),               0, diff(N3,epslon),               0, diff(N4,epslon),               0, diff(N5,epslon),               0, diff(N6,epslon),               0, diff(N7,epslon),               0, diff(N8,epslon),              0],
                     [              0,     diff(N1,ni),               0,     diff(N2,ni),               0,     diff(N3,ni),               0,     diff(N4,ni),               0,     diff(N5,ni),               0,     diff(N6,ni),               0,     diff(N7,ni),               0,    diff(N8,ni)],
                     [    diff(N1,ni), diff(N1,epslon),     diff(N2,ni), diff(N2,epslon),     diff(N3,ni), diff(N3,epslon),     diff(N4,ni), diff(N4,epslon),     diff(N5,ni), diff(N5,epslon),     diff(N6,ni), diff(N6,epslon),     diff(N7,ni), diff(N7,epslon),     diff(N8,ni), diff(N8,epslon)]
        
        ])

        D = self.D
        
        return t * (0.25) * (x5-x7) * (y6-y8) * integrate( B.T * D * B, (epslon,-1.0,1.0) , (ni,-1.0,1.0) )

    
    def calcular_kge(self,pontos,Kge):
        
        kl = self.kl
        
        ref_pontos = self.definir_ref_pontos(pontos)
        self.ref_pontos = ref_pontos
        
        #Introduzir a KL em Kge
        for i in range(len(ref_pontos)*ngl):
               
            for j in range(len(ref_pontos)*ngl):                   
                
                j1 = (i-i%ngl)/ngl
                k1 = (j-j%ngl)/ngl
                
                x = ref_pontos[j1]+i%ngl
                y = ref_pontos[k1]+j%ngl
                
                Kge[x,y] = Kge[x,y] + kl[i,j]                    

        return Kge
