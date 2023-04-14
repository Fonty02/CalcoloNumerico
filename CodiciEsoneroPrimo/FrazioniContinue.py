def calcolaValore(list: list,pos):

        if pos==len(list)-1 : return list[pos]
        else: return list[pos] + 1/calcolaValore(list,pos+1)



l=[3,4,12,4]
print(calcolaValore(l,0))