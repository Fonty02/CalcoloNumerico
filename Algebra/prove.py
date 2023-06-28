from numpy import *
from scipy import *

def Vandermon(x):
    righe,colonne=shape(x)
    V=zeros((righe,colonne))
    for i in range (righe):
        for j in range(colonne):
            V[i,j]=x[i]**j
    return V

def lagrange(x,y,z):
    yy=0
    n=shape(x)[0]
    for k in range (n):
        Lk=1
        for i in range(n):
            if i!=k:
                Lk*=(z-x[i])/(x[k]-x[i])
        yy+=Lk*y[k]
    return yy

def sostituzione_indietro(A,b):
    m,n=shape(A)
    x=zeros((n,1))
    tol=1e-15
    for i in range(m - 1, -1, -1):
        if abs(A[i, i]) < tol:
            raise ValueError('Matrice singolare')
        somma = 0
        for j in range(i + 1, n):
            somma = somma + A[i, j] * x[j]
        x[i] = (b[i] - somma) / A[i, i]
    return x

def sostituzione_avanti(A,b):
    m,n=shape(A)
    x=zeros((n,1))
    tol=1e-15
    for i in range(n):
        if abs(A[i,i])