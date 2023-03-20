from ZeriDIFunzione.zeri import *
import matplotlib.pyplot
import numpy



def showFunction(f,s_range = -3, e_range = 3):

    xlist = numpy.linspace(-7.5, 7.5, num=10000)
    ylist = f(xlist)
    fig = matplotlib.pyplot.figure(num=500, dpi=800)
    ax = fig.add_subplot(1, 1, 1)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('zero')
    ax.xaxis.set_label_coords(0, -0.05)
    ax.yaxis.set_label_coords(.001, 0)
    matplotlib.pyplot.plot(xlist, ylist,color="green")
    matplotlib.pyplot.title("f(x) = x - cos(x)",color="red",size=15,weight='bold')
    matplotlib.pyplot.ylabel("ASSE DELLE Y",color="red",size=15.,weight='bold')
    matplotlib.pyplot.xlabel("ASSE DELLE X",color="red",size=15,weight='bold')
    matplotlib.pyplot.grid()
    matplotlib.pyplot.show()


print("f(x) = x - cos(x)")
print("BISEZIONI ASSOLUTO = " + str(bisezioni_assoluto(f, 0.0000000001, pi / 4)))
print("BISEZIONI RELATIVO = " + str(bisezioni_relatvivo(f, 0.0000000001, pi / 4)))
print("BISEZIONI MISTO = " + str(bisezioni_misto(f, 0.0000000001, pi / 4)))
print("NEWTON = "+str(newton(f,pi/4)))
print("DIREZIONE COSTANTE ="+str(direzioneCostante(f,pi/4,0.5)))
print("NEWTON MODIFICATO ="+str(newtonModificato(f,pi/4)))
print("METODO DELLE SECANTI RELATIVO="+str(metodoDelleSecantiRelativo(f,8,4)))
print("METODO DELLE SECANTI MISTO ="+str(metodoDelleSecantiRelativo(f,8,4)))
showFunction(f)
