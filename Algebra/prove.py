import numpy
from numpy import *

def fattLu(A:numpy.array):
    n,m=numpy.shape(A)
    if n!=m: raise Exception("Matrice non quadrata")
    A=copy(A) #creo la copia locale dato che passa per riferimento
    L=numpy.ones((n,n))
    tol=1e-15
    for k in range(n-1):
        if abs(A[k,k])<tol: raise Exception("ELEMENTO PIVOTALE NULLO")
        for i in range (k+1,n):
            mik=-A[i,k]/A[k,k]
            L[i,k]=-mik
            for j in range (k+1,n):
                A[i,j]=A[i,j]+mik*A[k,j]
    U=triu(A)
    return L,U


def triangolareSUperiore(A:numpy.array,b:numpy.array):  #algoritmo di sost. all'indietro
    n=shape(b)[0]
    if shape(A)!=(n,n): raise Exception("NON VALIDO")
    x=numpy.zeros((shape(b)))
    for i in range (n): #parte dall'ultimo elemento (n-1) e deve arrivare al primo (quindi diminusce fino a -1 escluso, quindi 0)
        if abs(A[i,i])<1e-15: raise Exception("Matrice singolare")
        sum=0
        for j in range(i):
            sum+=A[i,j]*x[j]
        x[i]=(b[i]-sum)/A[i,i]
    return x



