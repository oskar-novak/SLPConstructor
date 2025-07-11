def logOps(Normalizer):
    n_qubits=Normalizer.shape[1]//2
    already_visited=[]
    V = np.block([
            [np.zeros((n_qubits, n_qubits)), np.eye(n_qubits)],
            [np.eye(n_qubits), np.zeros((n_qubits, n_qubits))]])
    logicalOperators = np.empty((0, Normalizer.shape[1]))
    for i in range(Normalizer.shape[0]):
        check=(Normalizer@V@Normalizer[i,:].T)%2
        if np.all(check == 0) or i in already_visited:
            continue
        indices_of_one = np.argwhere(check == 1)
        indices_of_one=indices_of_one.reshape(-1)
        already_visited.append(indices_of_one[0])
        
        logop1=Normalizer[i,:]
        logop2=Normalizer[indices_of_one[0],:]
        logop1=logop1[np.newaxis, :]
        logop2=logop2[np.newaxis, :]
        
        logicalOperators= np.append(logicalOperators, logop1, axis=0)
        logicalOperators= np.append(logicalOperators, logop2, axis=0)
        for j in range(i+1,Normalizer.shape[0]):
            if j==indices_of_one[0]:
                continue
            Normalizer[j,:]=(Normalizer[indices_of_one[0],:]*(Normalizer[i,:]@np.fft.fftshift(Normalizer[j,:]).transpose()%2)+ 
            Normalizer[j,:]) %2
            Normalizer[j,:]=((Normalizer[i,:]*(Normalizer[indices_of_one[0],:]@np.fft.fftshift(Normalizer[j,:]).transpose()%2)+
            Normalizer[j,:])%2)
           

    return np.unique(logicalOperators, axis=0)
