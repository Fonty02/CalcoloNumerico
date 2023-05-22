import numpy
from numpy import *

def sostituzioneIndietro(a,b):
    r,c=shape(a)
    if r!=c: raise ("Matrice non quadrata")
    x=zeros(r)
    tol=1e-15
    for i in range(r-1,-1,-1):
        if abs(a[i,i])<tol: raise ("Matrice singolare")
        sum=0
        for j in range (i+1,r):
            sum+=a[i,j]*x[j]
        x[i]=(b[i]-sum)/a[i,i]
    return x



def sostituzioneAvanti(a,b):
    r,c=shape(a)
    if r!=c: raise ("Matrice non quadrata")
    x=zeros(r)
    tol=1e-15
    for i in range(r):
        if abs(a[i,i])<tol: raise ("Matrice singolare")
        sum=0
        for j in range (i):
            sum+=a[i,j]*x[j]
        x[i]=(b[i]-sum)/a[i,i]
    return x

def fattLU(a):
    m,n=shape(a)
    if m!=n: raise("Matrice non quadrata")
    tol=1e-15
    a=copy(a)
    l=identity(n)
    for k in range(n-1):
        if abs(a[k,k])<tol: raise ("Matrice singolare")
        for i in range(k+1,n):
            mik=-a[i,k]/a[k,k]
            a[i]=a[i]+mik*a[k]
            l[i,k]=-mik
    return l,triu(a)

def eliminazioneGauss(a,b):
    m,n=shape(a)
    if m!=n: raise("Matrice non quadrata")
    tol=1e-15
    a=copy(a)
    b=copy(b)
    for k in range(n-1):
        if abs(a[k,k])<tol: raise("Matrice singolare")
        for i in range(k+1,n):
            mik=-a[i,k]/a[k,k]
            a[i]=a[i]+mik*a[k]
            b[i]=b[i]+mik*b[k]
    return sostituzioneIndietro(triu(a),b)

def massimoPivotParziale(a,b):
    m,n=shape(a)
    if m!=n: raise("Matrice non quadrata")
    a=copy(a)
    b=copy(b)
    for k in range(n-1):
        pivot=abs(a[k,k])
        s=k
        for i in range(k+1,n):
            if abs(a[i,k])>pivot:
                s=i
                pivot=abs(a[i,k])
        if s!=k:
            a[[k,s]]=a[[s,k]]
            b[s],b[k]=b[k],b[s]
        for i in range(k+1,n):
            mik=-a[i,k]/a[k,k]
            a[i]=a[i]+mik*a[k]
            b[i]=b[i]+mik*b[k]
    return sostituzioneIndietro(triu(a),b)



def fattLUconPivot(a):
    n,m=shape(a)
    if n!=m: raise("Matrice non quadrata")
    a=copy(a)
    l=identity(n)
    p=identity(n)
    tol=1e-15
    for k in range(n-1):
        if abs(a[k,k])<tol:
            i=k+1
            trovato=False
            while (i<n) and not trovato:
                if abs(a[i,k])>tol:
                    trovato=True
                    a[[i,k]]=a[[k,i]]
                    l[[i, k]] = l[[k, i]]
                    p[[i, k]] = p[[k, i]]
                else: i+=1
            if not trovato: raise("IMPOSSIBILE")
        for i in range(k+1,n):
            mik=-a[i,k]/a[k,k]
            a[i]=a[i]+mik*a[k]
            l[i,k]=-mik
    return p,l,triu(a)


def inversaLU(a):
    n,m=shape(a)
    if n!=m: raise("Matrice non quadrata")
    a=copy(a)
    i=identity(n)
    inv=identity(n)
    tol=1e-15
    for k in range(n-1):
        if abs(a[k,k])<tol: raise("Matrice singolare")
        for i in range (k+1,n):
            mik=-a[i,k]/a[k,k]
            a[i]=a[i]+mik*a[k]
            i[i]=i[i]+mik*i[k]
    u=triu(a)
    for i in range(n):
        inv[:,i]=sostituzioneIndietro(u,i[:,i])
    return inv

def riduzioneScalini(a):
    righe,colonne=shape(a)
    tol=1e-15
    a=copy(a)
    i=0
    j=0
    trovato=True
    while i<(righe-1) and j<(colonne) and trovato:
        trovato=abs(a[i,j])>tol
        while j<colonne and not trovato:
            h=i
            while h<righe and not trovato:
                if abs(a[h,j])>tol:
                    a[[h,i]]=a[[i,h]]
                    trovato=True
                else: h+=1
            if not trovato: j+=1
        for k in range(i+1,righe):
            mik=-a[k,j]/a[i,j]
            a[k]=a[k]+mik*a[i]
        i+=1
        j+=1
    return a

def rank(a):
    a=copy(a)
    a=riduzioneScalini(a)
    r,c=shape(a)
    count=0
    for i in range(r):
        if not (numpy.all((a[i]==0))):
            count+=1
    return count




A=array([[-3,-1,1,1],[-9,-1,-4,3],[-6,4,0,2]])
print(riduzioneScalini(A))
print(rank(A))
