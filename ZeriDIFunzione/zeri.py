from numpy import *


def bisezioni_assoluto(f, a, b, tol):
    """
    Metodo delle successive bisezioni
    ---------------------------------
    Sintassi
    -------------
        alpha=bisezioni(f,a,b,tol,itmax)

    Dati di input
    -----------------------
    f: funzione di cui calcolare uno zero
    [a,b]: intervallo di lavoro
    tol:  precisione richiesta

    -------
    Dati di input

    c: uno zero di f
    it: numero di iterate eseguite
    """
    itmax = ceil(log2((b - a) / tol))
    fa = f(a)
    fb = f(b)
    try:
        if fa * fb >= 0:
            raise Exception("Intervallo non valido (la funzione non cambia di segno)")
    except Exception:
        raise
    it = 0
    while it < itmax:
        it += 1
        c = (a + b) / 2
        fc = f(c)
        if fc == 0:
            return c, it
        elif fa * fc < 0:
            b = c
        else:
            a = c
            fa = fc
    if (b - a) > tol:
        print("Precisione non raggiunta, errore assoluto")
    return c, it


def bisezioni_relatvivo(f, a, b, tol, itmax):
    """
    Metodo delle successive bisezioni
    ---------------------------------
    Sintassi
    -------------
        alpha=bisezioni(f,a,b,tol,itmax)

    Dati di input
    -----------------------
    f: funzione di cui calcolare uno zero
    [a,b]: intervallo di lavoro
    tol:  precisione richiesta
    itmax: numero di iterate massimo

    -------
    Dati di input

    c: uno zero di f
    it: numero di iterate eseguite


    ----
    AVVERTENZE

    Non utilizzabile se a=0 OR b=0
    """
    fa = f(a)
    fb = f(b)
    try:
        if fa * fb >= 0:
            raise Exception("Intervallo non valido (la funzione non cambia di segno)")
    except Exception:
        raise
    it = 0
    while ((b - a) / min(abs(a), abs(b))) > tol and it < itmax:
        it += 1
        c = (a + b) / 2
        fc = f(c)
        if fc == 0:
            return c, it
        elif fa * fc < 0:
            b = c
        else:
            a = c
            fa = fc
    if (b - a) > tol:
        print("Precisione non raggiunta, errore relativo")
    return c, it


def bisezioni_misto(f, a, b, tol, itmax):
    """
    Metodo delle successive bisezioni
    ---------------------------------
    Sintassi
    -------------
        alpha=bisezioni(f,a,b,tol,itmax)

    Dati di input
    -----------------------
    f: funzione di cui calcolare uno zero
    [a,b]: intervallo di lavoro
    tol:  precisione richiesta
    itmax: numero di iterate massimo

    -------
    Dati di input

    c: uno zero di f
    it: numero di iterate eseguite
    """
    fa = f(a)
    fb = f(b)
    try:
        if fa * fb >= 0:
            raise Exception("Intervallo non valido (la funzione non cambia di segno)")
    except Exception:
        raise
    it = 0
    while ((b - a) / (1 + min(abs(a), abs(b)))) > tol and it < itmax:
        it += 1
        c = (a + b) / 2
        fc = f(c)
        if fc == 0:
            return c, it
        elif fa * fc < 0:
            b = c
        else:
            a = c
            fa = fc
    if (b - a) > tol:
        print("Precisione non raggiunta, errore misto")
    return c, it


def newton(f, x0, tol=1e-10, itmax=100):
    """
    Metodo delle successive bisezioni
    ---------------------------------
    Sintassi
    --------
      alpha=newton(f,x0,tol,itmax)

    Dati di input
    -------------
      f: funzione di cui calcolare uno zero

      x0: stima iniziale dello zero di f

      tol: precisione richiesta

      itmax: numero massimo di iterate consentite

    Dati di output
    --------------
      x1: approssimazione di uno zero di f(x)

      it: numero di iterate eseguite

    Autore: F. Iavernaro
    """
    arresto = False
    it = 0  # contatore di iterate
    while not (arresto) and it < itmax:
        it = it + 1
        x1 = x0 - f(x0) / f(x0, 1)  # 1 -f(x0)/f'(x0)
        arresto = abs(x1 - x0) < tol
        x0 = x1
    if not (arresto):
        print('Attenzione: precisione non raggiunta')
    return x1, it

def direzioneCostante(f,x0,m,tol=1e-10,itmax=1000):
    arresto=False
    it=0
    while not(arresto) and it<itmax:
        it+=1
        x1=x0-m*f(x0)
        arresto=abs(x1-x0)<tol
        x0=x1
    if not(arresto):
        print("Precisione non raggiunta")
    return x1,it
def f(x,ord=0):
    if ord==0:
        y=x-cos(x)
    elif ord==1:
        y=1+sin(x)
    else:
        print('ordine di derivazione non definito')
        return
    return y



