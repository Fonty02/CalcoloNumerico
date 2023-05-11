from matplotlib import pyplot
from scipy import *
from numpy import *


def Vandermonde(x):
    r=c=shape(x)[0]
    V=zeros((r,c))
    V[:,0]=ones((1,r))
    for i in range(r):
        for j in range(c):
            V[i,j]=x[i]**(j)
    return V

def solveVandermond(x,y):
    return linalg.solve(Vandermonde(x),y)


def lagrange(x,y,xx):
    """
    INPUT
    x: vettore dei nodi
    y: vettore delle ordinate
    xx: vettore di ascisse in cui valutare pn(x)

    OUTPUT
    yy: vettore delle valutazioni di pn(x) nelle ascisse xx(i)
    """
    n=shape(x)[0]  #numero di nodi
    yy=0
    for k in range(n):
        Lk=1
        for i in range(n):
            if i!=k: Lk*=((xx-x[i])/(x[k]-x[i]))  #xx è un vettore. Quindi a ogni elemento di xx sottraggo x[i] e ottengo
                                                  # un vettore. Poi CON *
                                                  # STO FACENDO IL PRODOTTO TRA GLI ELEMENTI DELLO STESSO INDICE
        yy+=(y[k]*Lk)
    return yy

def interpola(f,int,n):
    """
    INPUT
    f : funzione
    int: intervallo di interpolazione
    n: numero di nodi equidistanti che ricoprono l'intervallo int
    """
    x=linspace(int[0],int[1],n) #vettore dei nodi
    y=f(x) #vettore delle ordinate
    xx=linspace(int[0],int[1],100) #vettore di ascisse in cui valutare pn(x)
    yy=lagrange(x,y,xx) #vettore delle valutazioni di pn(x) nelle ascisse xx(i)
    fxx=f(xx) #valutazione della funzione f nelle ascisse xx(i)
    pyplot.plot(x,y,'bo',xx,fxx,'red',xx,yy,'green')
    pyplot.grid()
    pyplot.title="INTERPOLAZIONE DI LAGRANGE"
    pyplot.xlabel="ASSE DELLE X"
    pyplot.ylabel="ASSE DELLE Y"
    pyplot.show()

def f(x):
    return e**(-x)*sin(x)

print(interpola(f,[0,pi],100))