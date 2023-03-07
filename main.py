from ZeriDIFunzione.zeri import *
import matplotlib.pyplot
import numpy



def showFunction(f):
    xlist = numpy.linspace(-10, 10, num=10000)
    ylist = f(xlist)
    matplotlib.pyplot.figure(num=100, dpi=500)
    matplotlib.pyplot.plot(xlist, ylist, label="Funzione")
    matplotlib.pyplot.title("GRAFICO DI F")
    matplotlib.pyplot.ylabel("ASSE DELLE Y")
    matplotlib.pyplot.xlabel("ASSE DELLE X")
    matplotlib.pyplot.show()


print("ASSOLUTO = " + str(bisezioni_assoluto(f, 0.0000000001, pi / 2, 1e-10)))
print("RELATIVO = " + str(bisezioni_relatvivo(f, 0.0000000001, pi / 2, 1e-10, 10000)))
print("MISTO = " + str(bisezioni_misto(f, 0.0000000001, pi / 2, 1e-10, 10000)))

showFunction(f)
