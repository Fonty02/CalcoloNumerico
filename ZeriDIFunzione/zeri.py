from numpy import *


def f(x,ord=0):
    """

    :param x: valore dell'incognita da valutare
    :param ord: ordine della funzione
    :return: valore della funzione valutata in x
    """
    if ord==0:
        y=x-cos(x)
    elif ord==1:
        y=1+sin(x)
    else:
        print('ordine di derivazione non definito')
        return
    return y

def bisezioni_assoluto(f, a, b, tol=1e-10):
    """
    Metodo delle successive bisezioni con errore assoluto come criterio di stop
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


def bisezioni_relatvivo(f, a, b, tol=1e-10,itmax=100):
    """
    Metodo delle successive bisezioni con errore relativo come criterio di stop
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


def bisezioni_misto(f, a, b, tol=1e-10,itmax=100):
    """
    Metodo delle successive bisezioni con errore misto come criterio di stop
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
    Metodo di netwon per lo zero di funzioni
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
    """
    Metodo della direzione costante per lo zero di funzione
    ---------------------------------
    Sintassi
    --------
      alpha=direzioneCostante(f,x0,tol,itmax)

    Dati di input
    -------------
      f: funzione di cui calcolare uno zero

      x0: stima iniziale dello zero di f

      tol: precisione richiesta

      itmax: numero massimo di iterate consentite

      m: coefficiente angola della funzione iteratrice

    Dati di output
    --------------
      x1: approssimazione di uno zero di f(x)

      it: numero di iterate eseguite

    """
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



def newtonModificato(f, x0, tol=1e-10, itmax=100):
    """
    Metodo di Newton modificato (ovvero sfruttando la convergenza quadratica della direzione costante) per gli zeri di funzione
    ---------------------------------
    Sintassi
    --------
      alpha=direzioneCostante(f,x0,tol,itmax)

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

    """
    arresto=False
    it=0
    m=1/f(x0,1)
    while not(arresto) and it<itmax:
        it+=1
        x1=x0-m*f(x0)
        arresto=abs(x1-x0)<tol
        x0=x1
    if not(arresto):
        print("Precisione non raggiunta")
    return x1,it


def metodoDelleSecantiRelativo(f,x0,x1,tol=1e-10,itmax=100):
    """
    Metodo delle secanti per gli zeri di funzione
    ---------------------------------
    Sintassi
    --------
      alpha=secanti(f,x0,x1,tol,itmax)

    Dati di input
    -------------
      f: funzione di cui calcolare uno zero

      x0: stima iniziale dello zero di f
      x1: seconda stima dello zero

      tol: precisione richiesta

      itmax: numero massimo di iterate consentite


    Dati di output
    --------------
      x2: approssimazione di uno zero di f(x)

      it: numero di iterate eseguite

    """
    arresto=False
    it=0
    fx0=f(x0)
    fx1=f(x1)
    while not(arresto) and it<itmax:
        it+=1
        x2=x1 - (fx1/((fx1-fx0)/(x1-x0)))
        arresto=abs((x2-x1)/(x2))<tol
        x0=x1
        x1=x2
        fx0=fx1
        fx1=f(x1)
    if not(arresto):
        print("Precisione non raggiunta")
    return x1,it


def metodoDelleSecantiMisto(f,x0,x1,tol=1e-10,itmax=100):
    """
    Metodo delle secanti per gli zeri di funzione
    ---------------------------------
    Sintassi
    --------
      alpha=secanti(f,x0,x1,tol,itmax)

    Dati di input
    -------------
      f: funzione di cui calcolare uno zero

      x0: stima iniziale dello zero di f
      x1: seconda stima dello zero

      tol: precisione richiesta

      itmax: numero massimo di iterate consentite


    Dati di output
    --------------
      x2: approssimazione di uno zero di f(x)

      it: numero di iterate eseguite

    """
    arresto=False
    it=0
    fx0 = f(x0)
    fx1 = f(x1)
    while not(arresto) and it<itmax:
        it+=1
        x2 = x1 - (fx1 / ((fx1 - fx0) / (x1 - x0)))
        arresto=abs((x2-x1))/(1+abs(x2))<tol
        x0=x1
        x1=x2
        fx0 = fx1
        fx1 = f(x1)
    if not(arresto):
        print("Precisione non raggiunta")
    return x1,it



def metodoDelleSecantiAssoluto(f,x0,x1,tol=1e-10,itmax=100):
    """
    Metodo delle secanti per gli zeri di funzione
    ---------------------------------
    Sintassi
    --------
      alpha=secanti(f,x0,x1,tol,itmax)

    Dati di input
    -------------
      f: funzione di cui calcolare uno zero

      x0: stima iniziale dello zero di f
      x1: seconda stima dello zero

      tol: precisione richiesta

      itmax: numero massimo di iterate consentite


    Dati di output
    --------------
      x2: approssimazione di uno zero di f(x)

      it: numero di iterate eseguite

    """
    arresto=False
    it=0
    fx0 = f(x0)
    fx1 = f(x1)
    while not(arresto) and it<itmax:
        it+=1
        x2 = x1 - (fx1 / ((fx1 - fx0) / (x1 - x0)))
        arresto=abs(x2-x1)<tol
        x0=x1
        x1=x2
        fx0 = fx1
        fx1 = f(x1)
    if not(arresto):
        print("Precisione non raggiunta")
    return x1,it

