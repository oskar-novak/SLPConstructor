
import galois
import numpy as np
import itertools
from itertools import permutations


GF = galois.GF(2)

#Define Polynomials 

#See Galois README for Instructions


p1=galois.Poly([0,1,0],field=GF)
p2=galois.Poly([1,0,0],field=GF)
p0=galois.Poly([1],field=GF)


#!!!!ALWAYS GF2 for field!!!!

#Define Base Matrix
H=np.array([[p2,p1,p0],[p1,p2,p1]])


Sx3,Sz3,Gx3,Gz3,Lx3,Lz3=SLP(H,3)

print(np.all(lift(H@conjugate(findGenerator(H,3),3),3)==0))
print(np.all(Sx3@Sz3.transpose()%2==0))


#Tanner's Code -> SLP L=31 [775,128,20]* H is slightly under rank

p1=galois.Poly([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],field=GF)
p2=galois.Poly([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],field=GF)
p4=galois.Poly([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],field=GF)
p5=galois.Poly([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],field=GF)
p7=galois.Poly([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],field=GF)
p8=galois.Poly([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],field=GF)
p9=galois.Poly([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],field=GF)
p10=galois.Poly([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],field=GF)
p14=galois.Poly([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],field=GF)
p16=galois.Poly([0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],field=GF)
p18=galois.Poly([0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],field=GF)
p19=galois.Poly([0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],field=GF)
p20=galois.Poly([0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],field=GF)
p25=galois.Poly([0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],field=GF)
p28=galois.Poly([0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],field=GF)

HT=np.array([[p1,p2,p4,p8,p16],[p5,p10,p20,p9,p18],[p25,p19,p7,p14,p28]])


SxT,SzT,GxT,GzT,LxT,LzT=SLP(H,31)



print(np.all(lift(HT@conjugate(findGenerator(HT,31),31),31)==0))
print(np.all(SxT@SzT.transpose()%2==0))

#Can scale to very high lift sizes and codes



