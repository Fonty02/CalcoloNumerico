from numpy import *


def bisezioni_assoluto(f,a,b,tol):
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
    itmax=ceil(log2((b-a)/tol))
    fa=f(a)
    fb=f(b)
    try:
        if fa*fb>=0:
            raise Exception("Intervallo non valido (la funzione non cambia di segno)")
    except Exception:
        raise
    it=0
    while b-a>tol and it<itmax:
        it+=1
        c=(a+b)/2
        fc=f(c)
        if fc==0:
            return c, it
        elif fa*fc<0:
            b=c
        else:
            a=c
            fa=fc
    if (b-a)>tol:
        print("Precisione non raggiunta, errore assoluto")
    return c,it






def bisezioni_relatvivo(f,a,b,tol,itmax):
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
    fa=f(a)
    fb=f(b)
    try:
        if fa*fb>=0:
            raise Exception("Intervallo non valido (la funzione non cambia di segno)")
    except Exception:
        raise
    it=0
    while ((b-a)/min(abs(a),abs(b))) > tol and it<itmax:
        it+=1
        c=(a+b)/2
        fc=f(c)
        if fc==0:
            return c,it
        elif fa*fc<0:
            b=c
        else:
            a=c
            fa=fc
    if (b-a)>tol:
        print("Precisione non raggiunta, errore relativo")
    return c,it


def bisezioni_misto(f,a,b,tol,itmax):
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
    fa=f(a)
    fb=f(b)
    try:
        if fa*fb>=0:
            raise Exception("Intervallo non valido (la funzione non cambia di segno)")
    except Exception:
        raise
    it=0
    while ((b-a)/(1+min(abs(a),abs(b))))>tol and it<itmax:
        it+=1
        c=(a+b)/2
        fc=f(c)
        if fc==0:
            return c,it
        elif fa*fc<0:
            b=c
        else:
            a=c
            fa=fc
    if (b-a)>tol:
        print("Precisione non raggiunta, errore misto")
    return c,it


def f(x):
    return x-math.cos(x)