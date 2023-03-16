
from numpy import *
from pylab import *  # modulo per le funzioni grafiche




def cubica(z, ord=0):
    if ord == 0:
        w = z ** 3 - 1
    elif ord == 1:
        w = 3 * z ** 2
    else:
        print("secondo argomento non definito")
    return w


def newton(f, x0, tol=1e-10, itmax=1000):
    """
    Metodo ddi Newton
    ---------------------------------
    Sintassi
    --------
      [x1, it]=newton(f,x0,tol,itmax)

    Dati di input
    -------------
      f: funzione di cui calcolare uno zero

      x0: stima iniziale dello zero di f

      tol: precisione richiesta

      itmax: numero assimo di iterate consentite

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
        x1 = x0 - f(x0) / f(x0, 1)
        arresto = abs(x1 - x0) < tol
        x0 = x1
    if not (arresto):
        print('Attenzione: precisione non raggiunta')
    return x1, it



def frattale_cubica(xmin=-3.0, xmax=3.0, ymin=-3.0, ymax=3.0):
    """
    Determina le regioni di stabilita'
    nel piano complesso
    del metodo di Newton applicato alla
    funzione cubica f(z)=z**3-1



   Commenti: per ottenere il primo grafico, dalla console lanciare

             >>> frattale_cubica()

             In seguito, fare uno zoom su un particolare del frattale
             e acquisire i nuovi assi mediante il comando axis,
             ad esempio

             >>> w=axis()

             la variabile w è formalmente una lista di 4 elementi:
             con riferimento alla figura, i primi due rappresentano
             i valori minimo e massimo dell'ascissa e gli ultimi due
             i valori minimo e massimo dell'ordinata.
             Quindi è possibile rieseguire la funzione passando
             questi nuovi valori, al fine di ottenere
             un'immagine dettagliata della nuova regione selezionata:

             >>> frattale_cubica(w[0],w[1],w[2],w[3])

             Si può poi ripetere questa operazione per scoprire nuovi
             dettagli e verificare così la proprietà di autosimilarità

             Autore: Felice Iavernaro
               data: 20/10/2021

    """
    res_x = 1200  # num. di pixel sull'asse x
    res_y = 1200  # num. di pixel sull'asse y
    zero = [1, -1 / 2 + sqrt(3) / 2j, -1 / 2 - sqrt(3) / 2j]
    x = linspace(xmin, xmax, res_x)
    y = linspace(ymin, ymax, res_y)
    Z = zeros(shape=(res_x, res_y))
    for i in range(0, res_x):
        if mod(i, 100) == 0: print(i)
        for j in range(0, res_y):
            [alpha, it] = newton(cubica, complex(x[i], y[j]), 1e-4, 1000)
            if abs(alpha - zero[0]) < 1e-2:
                Z[i, j] = 0
            elif abs(alpha - zero[1]) < 1e-2:
                Z[i, j] = 0.5
            elif abs(alpha - zero[2]) < 1e-2:
                Z[i, j] = 1
            else:
                Z[i, j] = NaN
    figure(1)
    clf()
    imshow(transpose(Z), origin='lower', cmap=plt.cm.jet, interpolation='nearest', extent=(xmin, xmax, ymin, ymax))
    # plot(real(zero),imag(zero),'o')
    xlabel("Re(c)")
    ylabel("Im(c)")
    show()
    return Z
