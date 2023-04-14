import math
def approssimaPi(tol,itmax):
    it = 0
    arresto=False
    pi=1
    numerator = 0.0
    while (not arresto and it<itmax):
        it+=1
        numerator = math.sqrt(2.0 + numerator)
        pi *= (numerator / 2.0)
        if abs(pi- math.pi) < tol: arresto = True

    pi = (1.0 / pi) * 2.0
    return pi


def approssimaPi2(tol,itmax):
    it=0
    arresto=False
    sum=0
    denominatore=1
    op=0
    while it<itmax and not arresto:
        if (op%2)==0: sum+=1/denominatore
        else : sum-=1/denominatore
        denominatore+=2
        it+=1
        op+=1
        if abs(4*sum-math.pi)<tol: arresto=True
    return 4*sum




print(approssimaPi2(1e-16,100))
