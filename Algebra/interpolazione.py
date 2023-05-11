from scipy import *
from numpy import *


def Vandermond(x):
    r=shape(x)[0]
    c=r+1
    V=zeros((r,c))
    V[:,0]=ones((1,r))
    for i in range(r):
        for j in range(1,c):
            V[i,j]=x[i]**(j)
    return V