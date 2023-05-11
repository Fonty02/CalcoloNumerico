from numpy import *



def fattLU(A):
    r,c=shape(A)
    if r!=c : raise ("MATRICE NON QUADRATA")
    tol=1e-15
    L=zeros((r,c))
    A=copy(A)
    for k in range (r):
        if abs(A[k,k])<tol: raise("MATRICE SINGOLARE")
        for i in range(k+1,c):
            mik=-A[i,k]/A[k,k]
            L[i,k]=-mik
            for j in range(k+1,r):
                A[i,j]=A[i,j]+mik*A[k,j]
    return L,triu(A)


def eliminazioneGaussSemplice(A,b):
    r,c=shape(A)
    if r!=c: raise("MATRICE NON QUADRATA")
    tol=1e-15
    A=copy(A)
    b=copy(b)
    for k in range(r):
        if abs(A[k,k])<tol : raise ("MATRICE SINGOLARE")
        for i in range(k+1,c):
            mik=-A[i,k]/A[k,k]
            b[i]=b[i]+mik*b[k]
            for j in range(k+1,c):
                A[i,j]=A[i,j]+mik*A[k,j]
    return linalg.solve(triu(A),b)


def massimoPivotParziale(A,b):
    r, c = shape(A)
    if r != c: raise ("MATRICE NON QUADRATA")
    A = copy(A)
    b = copy(b)
    for k in range(r):
        pivot=abs(A[k,k])
        s=k
        for i in range(k+1,r):
            if abs(A[i,k])>pivot:
                s=i
                pivot=abs(A[i,k])
        if s!=k:
            A[[k,s]]=A[[s,k]]
            b[i],b[k]=b[k],b[i]
        for i in range(k + 1, c):
            mik = -A[i, k] / A[k, k]
            b[i] = b[i] + mik * b[k]
            for j in range(k + 1, c):
                A[i, j] = A[i, j] + mik * A[k, j]
    return linalg.solve(triu(A),b)



def fattLUconPivot(A):
    r,c=shape(A)
    if r!=c: raise ("MATRICE QUADRATA")
    L=zeros((r,c))
    P=identity(r)
    tol=1e-15
    A=copy(A)
    for k in range(r):
        if abs(A[k,k])<tol:
            trovato=False
            i=k+1
            while(i<r and not trovato):
                if abs(A[i,k])>tol:
                    A[[i,k]]=A[[k,i]]
                    L[[i,k]]=A[[k,i]]
                    P[[i,k]]=P[[k,i]]
                    trovato=True
                    break
                else: i+=1
            if not trovato: raise ("IMPOSSIBILE")
        for i in range(k + 1, c):
            mik = -A[i, k] / A[k, k]
            L[i, k] = -mik
            for j in range(k + 1, r):
                A[i, j] = A[i, j] + mik * A[k, j]
        return P, L, triu(A)




