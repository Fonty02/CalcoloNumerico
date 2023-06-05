from numpy import *
from scipy import *

def sostituzione_indietro(A,b):
    m,n=shape(A)
    if m!=n:
        print("Errore")
        return
    x=zeros((n,1))
    tol=1e-15
    for i in range(n-1,-1,-1):
        if abs(A[i,i])<tol:
            print("SIngolare")
            return
        sum=0
        for j in range(i+1,n):
            sum+=A[i,j]*x[j]
        x[i]=(b[i]-sum)/A[i,i]
    return x

def sostituzione_avanti(A,b):
    m,n=shape(A)
    if m!=n:
        print("Errore")
        return
    x=zeros((n,1))
    tol=1e-15
    for i in range(n):
        if abs(A[i,i])<tol:
            print("SIngolare")
            return
        sum=0
        for j in range(i):
            sum+=A[i,j]*x[j]
        x[i]=(b[i]-sum)/A[i,i]
    return x

def fattLU(A):
    m,n=shape(A)
    if m!=n:
        print("Matrice non quadrata")
        return
    tol=1e-15
    A=copy(A)
    L=identity(n)
    for k in range(n-1):
        if abs(A[k,k])<tol:
            print("Matrice singolare")
            return
        for i in range(k+1,n):
            mik=-A[i,k]/A[k,k]
            L[i,k]=-mik
            for j in range(k+1,n):
                A[i,j]=A[i,j]+mik*A[k,j]
    return L,triu(A)

def eliminazioneGauss(A,b):
    m,n=shape(A)
    if m!=n:
        print("Matrice non quadrata")
        return
    tol=1e-15
    A=copy(A)
    b=copy(b)
    for k in range(n-1):
        if abs(A[k,k])<tol:
            print("Matrice singolare")
            return
        for i in range(k+1,n):
            mik=-A[i,k]/A[k,k]
            b[i]=b[i]+mik*b[k]
            for j in range(k+1,n):
                A[i,j]=A[i,j]+mik*A[k,j]
    return linalg.solve(triu(A),b)

def massimoPivotParziale(A,b):
    m,n=shape(A)
    if m!=n:
        print("Matrice non quadrata")
        return
    A=copy(A)
    b=copy(b)
    for k in range(n-1):
        s=k
        pivot=abs(A[k,k])
        for i in range(k+1,n):
            tmp=abs(A[i,k])
            if tmp>pivot:
                s=i
                pivot=tmp
        if s!=k:
            A[[k,s]]=A[[s,k]]
            b[k],b[s]=b[s],b[k]
        for i in range(k+1,n):
            mik=-A[i,k]/A[k,k]
            b[i]=b[i]+mik*b[k]
            for j in range(k+1,n):
                A[i,j]=A[i,j]+mik*A[k,j]
    return linalg.solve(triu(A),b)

def fattLUconPivot(A):
    m,n=shape(A)
    if n!=m:
        raise "Non quadrata"
    A=copy(A)
    L=identity(n)
    P=identity(n)
    tol=1e-15
    for k in range(n-1):
        if abs(A[k,k])<tol:
            trovato=False
            for i in range(k+1,n):
                if abs(A[i,k])>tol:
                    A[[i,k]]=A[[k,i]]
                    L[[i,k]]=L[[k,i]]
                    P[[i,k]]=P[[k,i]]
                    trovato=True
                    break
            if not trovato:
                raise "Errore"
        for i in range(k+1,n):
            mik=-A[i,k]/A[k,k]
            L[i,k]=-mik
            for j in range(k+1,n):
                A[i,j]=A[i,j]+mik*A[k,j]
    return P,L,triu(A)

def inversaLU(A):
    m,n=shape(A)
    if m!=n:
        print("Matrice non quadrata")
        return
    tol=1e-15
    A=copy(A)
    I=identity(n)
    X=zeros((n,1))
    for k in range(n-1):
        if abs(A[k,k])<tol:
            print("Matrice singolare")
            return
        for i in range(k+1,n):
            mik=-A[i,k]/A[k,k]
            for j in range(k+1,n):
                A[i,j]=A[i,j]+mik*A[k,j]
            for j in range(k , n):
                I[i, j] = I[i, j] + mik * I[k, j]
    U=triu(A)
    for i in range(n):
        X[:,i]=linalg.solve(U,I[:,i])
    return X

def riduzioneScalini(A):
    righe,colonne=shape(A)
    A=copy(A)
    tol=1e-15
    trovato=True
    i=0
    j=0
    while i<righe-1 and j<colonne and trovato:
        trovato=abs(A[i,j])>tol
        while j<colonne and not trovato:
            h=i
            while h<righe and not trovato:
                if abs(A[h,j])>tol:
                    trovato=True
                    A[[h,i]]=A[[i,h]]
                    break
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
    for i in range(shape(C)[0]):
        if sum(C[i])!=0: count+=1
    return count

def norma(A,c):
    A=copy(A)
    if c==1: A=transpose(A)
    s=sum(abs(A[0]))
    for i in range(1,shape(A)[0]):
        if sum(abs(A[i]))>s:
            s=sum(abs(A[i]))
    return s

def potenze(A,y0,tol=1e-15,kmax=500):
    z0=y0/linalg.nomr(y0)
    sigma0=0
    k=0
    arresto=False
    while k<kmax and not arresto:
        t=dot(A,z0)
        z1=t/linalg.norm(t)
        sigma1=sum(t*z0)
        Er=abs(sigma1 - sigma0)/abs(sigma1)
        arresto= abs(Er)<tol
        sigma0=sigma1
        z0=z1
    if not arresto:
        raise "Error"
    return sigma1,z1


def Vandermonde(x):
    righe,colonne=shape(x)
    V=zeros((righe,righe))
    for i in range(righe):
        for j in range(righe):
            V[i,j]=x[i]**j
    return V


def lagrande(x,y,xx):
    yy=0
    n=shape(x)[0]
    for k in range(n):
        Lk=1
        for i in range(n):
            if i!=k: Lk*=(xx-x[i])/(x[k]-x[i])
        yy+=Lk*y[k]
    return yy