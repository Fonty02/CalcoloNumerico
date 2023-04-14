def collatz(nInput):
    l=[]
    print(nInput)
    while (nInput!=1):
        l.append(nInput)
        if nInput % 2 == 0:
            nInput=(nInput // 2)        # // serve per prendere la parte intera perchè Python quando vede il / trasforma subito in reale
        else : nInput=(3 * nInput + 1)
    l.append(nInput)
    return l

nInput = int(input())
print(collatz(nInput))