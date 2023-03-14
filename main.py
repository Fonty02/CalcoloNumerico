from ZeriDIFunzione.zeri import *
import matplotlib.pyplot
import numpy



def showFunction(f,s_range = -3, e_range = 3):
    plot = numpy.linspace(s_range, e_range)

    # setting the axes at the centre
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # plot the function
    matplotlib.pyplot.plot(plot, f(plot),'r')
    matplotlib.pyplot.grid()
    # show the plot
    matplotlib.pyplot.show()


print("ASSOLUTO = " + str(bisezioni_assoluto(f, 0.0000000001, pi / 2, 1e-10)))
print("RELATIVO = " + str(bisezioni_relatvivo(f, 0.0000000001, pi / 2, 1e-10, 10000)))
print("MISTO = " + str(bisezioni_misto(f, 0.0000000001, pi / 2, 1e-10, 10000)))
print("NEWTON = "+str(newton(f,pi/2,1e-10, 10000)))
print("DIREZIONE COSTANTE ="+str(direzioneCostante(f,pi/2,0.5)))

#showFunction(f)
