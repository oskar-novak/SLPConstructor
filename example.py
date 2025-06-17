import galois
import numpy as np
GF = galois.GF(2)

#Define Polynomials 

#See Galois README for Instructions
p1=galois.Poly([1,1,1],field=GF)
p2=galois.Poly([0,1,1],field=GF)
p3=galois.Poly([0,1,0],field=GF)
p4=galois.Poly([1,0,0],field=GF)
p5=galois.Poly([1],field=GF)
p6=galois.Poly([0],field=GF)

#!!!!ALWAYS GF2 for field!!!!

#Define Base Matrix
H=np.array([[p4,p3,p5],[p3,p4,p3]])


Hx,Hz,Gx,Gz,Lx,Lz=SLP(H,3)




