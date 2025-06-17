def SLP(H,Lift_Size):
    #find Generator Matrix

    G=findGenerator(H,Lift_Size)

    GF = galois.GF(2)

    p1=galois.Poly([1],field=GF)
    p0=galois.Poly([0],field=GF)

    #construct Code
    baseIdentity=np.empty((H.shape[1],H.shape[1]),dtype=object)

    for i in range(H.shape[1]):
        for j in range(H.shape[1]):
            if i==j:
                baseIdentity[i,j]=p1
            else:
                baseIdentity[i,j]=p0
    
    
    Sx=lift(np.kron(H,G),Lift_Size)
    Sz=lift(np.kron(G,H),Lift_Size)

    Gx=lift(np.kron(H,baseIdentity),Lift_Size)
    Gz=lift(np.kron(baseIdentity,H),Lift_Size)
    
    Lx=lift(np.kron(baseIdentity,G),Lift_Size)
    Lz=lift(np.kron(G,baseIdentity),Lift_Size)

    return Sx, Sz, Gx, Gz, Lx, Lz

    
