# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 08:16:36 2020

@author: felix
"""

from numpy import * 
from scipy import *

def laplace(A):
    """
    Regola di Laplace per il calcolo del determinante di una matrice quadrata.
    Sviluppo di Laplace lungo la prima riga
    """
    [m,n]=shape(A)
    if n==1:
        d=A[0,0]
    else:
        d=0
        for j in range(0,n):
            A1j=delete(A,0,axis=0)
            A1j=delete(A1j,j,axis=1)
            d=d+(-1)**(j)*A[0,j]*laplace(A1j)
    return d

def triang_sup(A,b):
    """
    Algoritmo di sostituzione all'indietro
    per la risoluzione dei sistemi lineari
    triangolari superiori
    
    --------------------------------------

    INPUT
    --------------------------------------
    A: matrice dei coefficienti triangolare superiore
    
    b: vettore dei termini noti
    
    OUTPUT
    --------------------------------------
    x: vettore soluzione del sistema Ax=b
    """
    [m,n]=shape(A)# oppure n=len(b)
    x=zeros(shape=(n,1))# preallochiamo la memoria per x
    tol=1e-15
    for i in range(n-1,-1,-1):
        if abs(A[i,i])<tol:
            raise ValueError('Matrice singolare')
        else:
            somma=0
            for j in range(i+1,n):
                somma=somma+A[i,j]*x[j]
            x[i]=(b[i]-somma)/A[i,i]
    return x


def fattlu(A):
    """
    Fattorizzazione LU di una matrice
    sintassi:
        [L,U]=fattlu(A)

    INPUT
    A: matrice da fattorizzare
    OUTPUT
    L: matrice triangolare inferiore speciale
    U: matrice triangolare superiore    
    """
    [m,n]=shape(A)
    if m!=n:
        print("matrice non quadrata")
    A=copy(A)# altrimenti sovrascriviamo la A nella shell
    tol=1e-15
    L=identity(n)
    for k in range(0,n-1):
        if abs(A[k,k])<tol:
            print("minore principale nullo")
            return
        for i in range(k+1,n):
            mik=-A[i,k]/A[k,k]
            for j in range(k+1,n):
                A[i,j]=A[i,j]+mik*A[k,j]
            L[i,k]=-mik
    U=triu(A) # estrae la parte triang. sup. di A
    return L,U             
    
    
    
    
    
    
    
    
    
    
    
    
    
    