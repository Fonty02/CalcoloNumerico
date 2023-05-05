from scipy import *
from numpy import *

def sommaMatriciale(A, B):
    """
    Somma matriciale
    ----------------
    INPUT
    A: matrice m x n
    B: matrice m x n
    ----------------
    OUTPUT
    C: matrice m x n
    """
    [m, n] = shape(A)
    [p, q] = shape(B)
    if m != p or n != q:
        raise ValueError('Dimensioni non compatibili')
    C = zeros(shape=(m, n), dtype=type(A[0, 0]))
    for i in range(0, m):
        for j in range(0, n):
            C[i, j] = A[i, j] + B[i, j]
    return C


def prodottoMatriciale(A, B):
    """
    Prodotto matriciale
    -------------------
    INPUT
    A: matrice m x n
    B: matrice n x p
    -------------------
    OUTPUT
    C: matrice m x p

    """
    [m, n] = shape(A)
    [p, q] = shape(B)
    if n != p:
        raise ValueError('Dimensioni non compatibili')
    C = zeros(shape=(m, q), dtype=type(A[0, 0]))
    for i in range(0, m):
        for j in range(0, q):
            for k in range(0, n):
                C[i, j] += A[i, k] * B[k, j]
    return C


def trasposta(A):
    """
    Trasposta di una matrice
    ------------------------
    INPUT
    A: matrice m x n
    ------------------------
    OUTPUT
    C: matrice n x m
    """
    [m, n] = shape(A)
    C = zeros(shape=(n, m), dtype=type(A[0, 0]))
    for i in range(0, m):
        for j in range(0, n):
            C[j, i] = A[i, j]
    return C


def sottoMatricePrincipaleDiTesta(A, k):
    """
    Sotto matrice quadrata di testa
    -------------------------------
    INPUT
    A: matrice m x n
    k: intero
    -------------------------------
    OUTPUT
    C: matrice k x k
    """
    [m, n] = shape(A)
    if k > m or k > n or k < 1:
        raise ValueError('Dimensioni non compatibili')
    C = zeros(shape=(k, k), dtype=type(A[0, 0]))
    for i in range(0, k):
        for j in range(0, k):
            C[i, j] = A[i, j]
    return C


def laplace(A):
    """
    Regola di Laplace per il calcolo del determinante di una matrice quadrata.
    Sviluppo di Laplace lungo la prima riga
    """
    [m, n] = shape(A)
    if n == 1:
        d = A[0, 0]
    else:
        d = 0
        for j in range(0, n):
            A1j = delete(A, 0, axis=0)
            A1j = delete(A1j, j, axis=1)
            d = d + (-1) ** (j) * A[0, j] * laplace(A1j)
    return d


def matriceAggiunta(A):
    """
    Calcolo della matrice aggiunta di una matrice quadrata
    """
    [m, n] = shape(A)
    if m != n:
        raise ValueError('Matrice non quadrata')
    B = zeros(shape=(m, n), dtype=type(A[0, 0]))
    for i in range(0, m):
        for j in range(0, n):
            Aij = delete(A, i, 0)
            Aij = delete(Aij, j, 1)
            Aji = transpose(Aij)
            B[j, i] = (-1) ** (i + j) * laplace(Aji)
    return B


def matriceInversa(A):
    """
    Calcolo della matrice inversa di una matrice quadrata
    """
    [m, n] = shape(A)
    if m != n:
        raise ValueError('Matrice non quadrata')
    d = laplace(A)
    B = matriceAggiunta(A)
    return B / d


def metodoCramer(A, b):
    """
    Metodo di Cramer per la risoluzione di sistemi lineari
    """
    [m, n] = shape(A)
    if m != n:
        raise ValueError('Matrice non quadrata')
    d = laplace(A)
    if abs(d) < 1e-15:
        raise ValueError('Matrice singolare')
    x = zeros(shape=(n, 1))
    for i in range(0, n):
        Ai = copy(A)  # copia di A
        Ai[:, i] = b  # sostituzione della i-esima colonna di A con b
        x[i] = laplace(Ai) / d  # calcolo della i-esima componente di x
    return x


def isTriangSup(A):
    """
    Verifica se una matrice è triangolare superiore
    """
    m, n = shape(A)
    if m != n:
        raise ValueError("LA MATRICE NON QUADRATA")
    for i in range(1, n):
        for j in range(0, i):
            if abs(A[i, j]) > 1e-15: return False
    return True


def triang_sup(A, b):
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
    if (not isTriangSup(A)):
        raise ValueError('Matrice non triangolare superiore')
    [m, n] = shape(A)  # oppure n=len(b)
    x = zeros(shape=(n, 1))  # preallochiamo la memoria per x
    tol = 1e-15
    for i in range(n - 1, -1, -1):
        if abs(A[i, i]) < tol:
            raise ValueError('Matrice singolare')
        somma = 0
        for j in range(i + 1, n):
            somma = somma + A[i, j] * x[j]
        x[i] = (b[i] - somma) / A[i, i]
    return x


def isTriangInf(A):
    """
    Verifica se una matrice è triangolare inferiore
    """
    m, n = shape(A)
    if m != n:
        raise ValueError("LA MATRICE NON QUADRATA")
    for i in range(0, n):
        for j in range(1, i):
            if abs(A[-i, -j]) > 1e-15: return False
    return True


def triang_inf(A, b):
    """
    Algoritmo di sostituzione in avanti
    per la risoluzione dei sistemi lineari
    triangolari inferiore

    --------------------------------------

    INPUT
    --------------------------------------
    A: matrice dei coefficienti triangolare inferiori

    b: vettore dei termini noti

    OUTPUT
    --------------------------------------
    x: vettore soluzione del sistema Ax=b
    """
    if (not isTriangInf(A)):
        raise ValueError('Matrice non triangolare inferiore')
    [m, n] = shape(A)  # oppure n=len(b)
    x = zeros(shape=(n, 1))  # preallochiamo la memoria per x
    tol = 1e-15
    for i in range(0, n):
        if abs(A[i, i]) < tol:
            raise ValueError('Matrice singolare')
        else:
            somma = 0
            for j in range(0, i):
                somma = somma + A[i, j] * x[j]
            x[i] = (b[i] - somma) / A[i, i]
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
    [m, n] = shape(A)
    if m != n:
        raise ValueError("matrice non quadrata")
    A = copy(A)  # altrimenti sovrascriviamo la A nella shell, la A a sinistra diventa variabile locale perchè si ha passaggio per riferimento
    tol = 1e-15
    L = identity(n)  #prendo la matrice identica di ordine n
    for k in range(0, n - 1):  #gli n-1 passi da fare (per calcolare A2, A3,...,An=U)
        if abs(A[k, k]) < tol:  #controllo che gli elementi della diagonale non siano nulli
            raise ValueError ("minore principale nullo")
        for i in range(k + 1, n):
            mik = -A[i, k] / A[k, k] #calcolo i moltiplicatori della matrice elementare di Gauss
            for j in range(k + 1, n):
                A[i, j] = A[i, j] + mik * A[k, j]  #nuova i-esima riga = vecchia i-esima + moltiplicatore per k-esami (dove sta l elemento pivotale)
            L[i, k] = -mik #creo gia la triangolare inferiore speciale
    U = triu(A)  # estrae la parte triang. sup. di A -> dato che non annullo a mano gli elementi delle colonne
    return L, U


def determinanteLU(A):
    L, U = fattlu(A)
    return laplace(U)


def risoluzioneSistemaLU(A, b):
    """
    Risoluzione di un sistema lineare sfruttando la fattorizzazioen LU
    """
    L, U = fattlu(A)
    y = triang_inf(L, b)
    x = triang_sup(U, y)
    sol = zeros(shape(x)[0])
    for i in range(0, shape(x)[0]):
        sol[i] = x[i, 0]
    return sol

def eliminazioneGauss(A,b):
    n,m=shape(A)
    if n!=m: raise Exception("Matrice non quadrata")
    A=copy(A)
    b=copy(b)
    tol=1e-15
    for k in range (n-1):
        if abs(A[k,k])<tol: raise Exception("Matrice singolare")
        for i in range(k+1,n):
            mik=-A[i,k]/A[k,k]
            b[i]=b[i]+mik*b[k]
            for j in range(k+1,n):
                A[i,j]=A[i,j]+mik*A[k,j]
    U=triu(A)
    return triang_sup(U,b)

def massimoPivotParziale(A,b):
    n, m = shape(A)
    if n != m: raise Exception("Matrice non quadrata")
    A = copy(A)
    b = copy(b)
    tol = 1e-15
    for k in range(n - 1):
        s=k
        pivot=abs(A[k,k])
        for i in range(k+1,n):
            if abs(A[i,k])>pivot:
                s=i
                pivot=abs(A[i,k])
        if s!=k:
            A[[k,s]]=A[[s,k]]
            b[k],b[s]=b[s],b[k]
        for i in range(k + 1, n):
            mik = -A[i, k] / A[k, k]
            b[i] = b[i] + mik * b[k]
            for j in range(k + 1, n):
                A[i, j] = A[i, j] + mik * A[k, j]
    U = triu(A)
    return triang_sup(U, b)

def fattLUconPivot(A):
    n, m = shape(A)
    if n != m: raise Exception("Matrice non quadrata")
    tol = 1e-15
    A = copy(A)
    L = identity(n)
    P = identity(n)
    for k in range(n - 1):
        if abs(A[k, k]) < tol:
            trovato = False
            i = k + 1
            while (not trovato and i < n):
                if abs(A[i, k]) > tol:
                    A[[i, k]] = A[[k, i]]
                    L[[i, k]] = L[[k, i]]
                    P[[i, k]] = P[[k, i]]
                    trovato = True
                    break
                else:
                    i += 1
            if not trovato: raise Exception("ELEMENTI SOTTOSTANTI AL PIVOT NULLI")
        for i in range(k + 1, n):
            mik = -A[i, k] / A[k, k]
            L[i, k] = -mik
            for j in range(k + 1, n):
                A[i, j] = A[i, j] + mik * A[k, j]
    U = triu(A)
    return P, L, U

def inversaLU(A):
    n,m=shape(A)
    if n!=m: raise Exception("Matrice non quadrata")
    A=copy(A)
    I=identity(n)
    inv=zeros((n,n))
    mik=0
    for k in range(0,n-1):
        if abs(A[k,k])<1e-15: raise Exception("Impossibile continuare")
        for i in range (k+1,n):
            mik=-A[i,k]/A[k,k]
            for j in range(k+1,n):
                A[i, j] = A[i, j] + mik * A[k, j]
            for j in range(n):
                I[i, j] = I[i, j] + mik * I[k, j]
    U=triu(A)
    for i in range(n):
        inv[:,i]=linalg.solve(U,I[:,i])
    return inv



def riduzioneScalini(A):
    A =copy(A)
    righe, colonne = shape(A)
    tol = 1e-15
    j = 0
    i = 0
    trovato=True
    while i<righe -1 and trovato:
        trovato= abs(A[i,j])>tol
        while j<colonne and not trovato:
            h = i
            while h<righe and not trovato:
                if abs(A[h,j])>tol:
                    trovato=True
                    A[[h,i]]=A[[i,h]]
                    break
                h+=1
            if not trovato: j+=1
        #calcoli
        if not trovato: raise Exception ("Impossibile proseguire")
        for k in range(i+1,righe):
            A[k]=A[k]-A[i]*(A[k,j]/A[i,j])
        i+=1
        j+=1
    return triu(A)

def rank(A):
    C=riduzioneScalini(A)
    count=0
    r,c=shape(C)
    vuoto=zeros((1,c))
    for i in range(r):
        vuota=True
        for j in range(c):
            if abs(C[i,j])>1e-15:
                vuota=False
                break
        if not vuota:
            count=count+1
    return count










A = array([[1,-3,-1,1], [3,-9,-1,-4], [-2,6,4,0]])
C=riduzioneScalini(A)
print(C)
print(rank(C))
exit(0)
