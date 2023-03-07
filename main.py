from ZeriDIFunzione.zeri import *
print("ASSOLUTO = "+str(bisezioni_assoluto(f,0.0000000001,pi/2,1e-17)))
print("RELATIVO = "+str(bisezioni_relatvivo(f,0.0000000001,pi/2,1e-17,10000)))
print("MISTO = "+str(bisezioni_misto(f,0.0000000001,pi/2,1e-17,10000)))



