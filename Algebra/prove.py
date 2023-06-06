from numpy import *
from scipy import *

def Vandermonde(x):
    m,n=shape(x)
    V=zeros((m,m))
    for i in range(m):
        for j in range(n):
            V[i,j]=x[i]**j
    return V

def lagrange(x,y,xx):
    yy=0
    n=shape(x)[0]
    for k in range(n):
        Lk=1
        for i in range(n):
            if i!=k: Lk*=(xx-x[i])/(x[k]-x[i])
        yy+=Lk*y[k]
    return yy

def sostituzione_indietro(A,b):
    n,m=shape(A)
    if n!=m:
        return
    tol=1e-15
    x=zeros((n,1))
    for i in range(n-1,-1,-1):
        if abs(A[i,i])<tol:
            return
        somma=0
        for j in range(i+1,n):
            somma+=A[i,j]*x[j]
        x[i]=(b[i]-somma)/A[i,i]
    return x

def sostituzione_avanti(A,b):
    n,m=shape(A)
    if n!=m:
        return
    tol=1e-15
    x=zeros((n,1))
    for i in range(n):
        if abs(A[i,i])<tol:
            return
        somma=0
        for j in range(i):
            somma+=A[i,j]*x[j]
    return x

def fattLU(A):
    n,m=shape(A)
    tol=1e-15
    A=copy(A)
    L=identity(n)
    for k in range(n-1):
        if abs(A[k,k])<tol:
            return
        for i in range(k+1,n):
            mik=-A[i,k]/A[k,k]
            L[i,k]=-mik
            for j in range(k+1,n):
                A[i,j]=A[i,j]+mik*A[k,j]
    return L,triu(A)

def eliminazioneGauss(A,b):
    n,m=shape(A)
    tol=1e-15
    A=copy(A)
    b=copy(b)
    for k in range(n-1):
        if abs(A[k,k])<tol:
            return
        for i in range (k+1,n):
            mik=-A[i,k]/A[k,k]
            b[i]=b[i]+mik*b[k]
            for j in range(k+1,n):
                A[i,j]=A[i,j]+mik*A[k,j]
    return sostituzione_indietro(triu(A),b)

def massimoPivotParziale(A,b):
    n,m=shape(A)
    A=copy(A)
    b=copy(b)
    tol=1e-15
    for k in range (n-1):
        pivot=abs(A[k,k])
        s=k
        for i in range(k+1,n):
            tmp=abs(A[i,k])
            if tmp>pivot:
                pivot=tmp
                s=i
        if s!=k:
            A[[k,s]]=A[[s,k]]
            b[k],b[s]=b[s],b[k]
    if abs(A[k,k])<tol:
        return
    for i in range(k+1,n):
        mik=-A[i,k]/A[k,k]
        b[i]=b[i]+mik*b[k]
        for j in range(k+1,n):
            A[i,j]=A[i,j]+mik*A[k,j]
    return sostituzione_indietro(triu(A),b)

def GaussTridiagonale(A,b):
    A=copy(A)
    tol=1e-15
    n,m=shape(A)
    b=copy(b)
    for k in range(n-1):
        if abs(A[k,k])<tol:
            return
        mik=-A[k+1,k]/A[k,k]
        b[k+1]=b[k+1]+mik*b[k]
        A[k,+1,k]=A[k+1,k]+mik*A[k,k+1]

def fattLUconPivot(A):
    n,m=shape(A)
    tol=1e-15
    A=copy(A)
    L=identity(n)
    P=copy(L)
    for k in range(n-1):
        trovato=False
        if abs(A[k,k])<tol:
            h=k+1
            while h<n and not trovato:
                if abs(A[h,k])>tol:
                    trovato=True
                    A[[k,h]]=A[[h,k]]
                    L[[k,h]]=A[[h,k]]
                    P[[k,h]]=P[[h,k]]
                h+=1
        if not trovato:
            return
        for i in range(k+1,n):
            mik=-A[i,k]/A[k,k]
            L[i,k]=-mik
            for j in range(k+1,n):
                A[i,j]=A[i,j]+mik*A[k,j]
    return P,L,triu(A)

def inversaLU(A):
    A=copy(A)
    n,m=shape(A)
    I=identity(n)
    X=zeros((n,n))
    tol=1e-15
    for k in range(n-1):
        if abs(A[k,k])<tol:
            return
        for i in range(k+1,n):
            mik=-A[i,k]/A[k,k]
            for j in range(k+1,n):
                A[i,j]=A[i,j]+mik*A[k,j]
            for j in range(k,n):
                I[i,j]=I[i,j]+mik*I[k,j]
    U=triu(A)
    for i in range(n):
        X[:,i]=sostituzione_indietro(U,I[:,i])
    return X

def riduzioneScalini(A):
    A=copy(A)
    righe,colonne=shape(A)
    i=0
    j=0
    trovato=True
    tol=1e-15
    while i<righe-1 and j<colonne and trovato:
        trovato=abs(A[i,j])>tol
        while j<colonne and not trovato:
            h=i
            while h<righe and not trovato:
                if abs(A[h,j])>tol:
                    trovato=True
                    A[[h,i]]=A[[i,h]]
                else: h+=1
            if not trovato: j+=1
        if not trovato: return triu(A)
        for k in range(i+1,righe):
            mik=-A[k,j]/A[i,j]
            A[k]=A[k]+mik*A[i]
        i+=1
        j+=1
    return triu(A)

def rank(A):
    C=riduzioneScalini(A)
    count=0
    righe=shape(C)[0]
    for i in range(righe):
        if sum(C[i])!=0:
            count+=1
    return count

def norma(A,c):
    if c=="1":
        A=transpose(A)
    righe,colonne=shape(A)
    massimo=sum(A[0])
    for i in range(1,righe):
        somma=sum(A[i])
        if somma>massimo:
            massimo=somma
    return massimo

def potenze(A,y0,tol=1e-10,kmax=500):
    arresto=False
    it=0
    z0=y0/linalg.norm(y0)
    sigma0=0
    while not arresto and it<kmax:
        t=dot(A,z0)
        z1=t/linalg.norm(t)
        sigma1=sum(t*z0)
        if (abs(sigma1-sigma0)/abs(sigma1))<tol:
            arresto=True
        it+=1
        z0=z1
        sigma0=sigma1
    if not arresto:
        return
    else: return sigma1,z1

