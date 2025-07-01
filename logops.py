def symplecticMatrix(x, y): #compact way to find commutativity relations between Pauli Ops
    z=np.vstack((x,y))
    num_rows, num_cols = z.shape
    n_qubits = num_cols // 2
    V = np.block([
            [np.zeros((n_qubits, n_qubits)), np.eye(n_qubits)],
            [np.eye(n_qubits), np.zeros((n_qubits, n_qubits))]
    ])
    return (x @ V @ y.T) % 2




def logOps(Normalizer):
    logicalOperators = np.empty((0, Normalizer.shape[1]))
    for i in range(Normalizer.shape[0]):
        check=symplecticMatrix(Normalizer,Normalizer[i,:])
        if np.all(check == 0):
            continue
        indices_of_one = np.argwhere(check == 1)
        indices_of_one=indices_of_one.reshape(-1)
        
        
        logop1=Normalizer[i,:]
        logop2=Normalizer[indices_of_one[0],:]
        logop1=logop1[np.newaxis, :]
        logop2=logop2[np.newaxis, :]
        
        logicalOperators= np.append(logicalOperators, logop1, axis=0)
        logicalOperators= np.append(logicalOperators, logop2, axis=0)
        for j in range(i,Normalizer.shape[0]):
            if j==i or j==indices_of_one[0]:
                continue
            Normalizer[j,:]=(Normalizer[indices_of_one[0],:]*symplecticMatrix(Normalizer[i,:],Normalizer[j,:])+Normalizer[j,:]) %2
            Normalizer[j,:]=((Normalizer[i,:]*symplecticMatrix(Normalizer[indices_of_one[0],:],Normalizer[j,:])+
            Normalizer[j,:])%2)
           

    return np.unique(logicalOperators, axis=0)
            
