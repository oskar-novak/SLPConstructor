def gf2elim(M):

    m,n = M.shape

    i=0
    j=0

    while i < m and j < n:
        
        # find value and index of largest element in remainder of column j
        k = np.argmax(M[i:, j]) +i

        # swap rows
      
        temp = np.copy(M[k])
        M[k] = M[i]
        M[i] = temp


        aijn = M[i, j:]

        col = np.copy(M[:, j]) #make a copy otherwise M will be directly affected

        col[i] = 0 #avoid xoring pivot row with itself

        flip = np.outer(col, aijn)

        M[:, j:] = M[:, j:] ^ flip

        i += 1
        j +=1

    return M

def rankMat(A):
    
    A=A.tolist()
    

    n=len(A[0])
    
    rank = 0
    for col in range(n):
        
        j=0
        
        rows = []
        
        while j<len(A):
            
            if A[j][col] == 1:
                
                rows += [j]
            j+=1
            
        if len(rows) >= 1:
            
            for c in range(1,len(rows)):
                
                for k in range(n):
                    
                    A[rows[c]][k] = (A[rows[c]][k] + A[rows[0]][k])%2
                    
            A.pop(rows[0])
            
            rank += 1
            
    for row in A:
        
        if sum(row)>0:
            
            rank += 1
            
    return rank  
