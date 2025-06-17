import galois
import numpy as np
import itertools
from itertools import permutations


def get_matrix_minor(arr, i, j):
    
    return np.delete(np.delete(arr, i, axis=0), j, axis=1)




def getMatrixDeterminant(m):
    
    GF=galois.GF(2)
    
    #base case for 2x2 matrix
    if len(m) == 2:
        
        return m[0,0]*m[1,1]-m[0,1]*m[1,0]

    determinant = galois.Poly([0],field=GF)
    
    for c in range(len(m)):
        
        minor=m[0,c]*getMatrixDeterminant(get_matrix_minor(m,0,c),GF.order)
        
        determinant= determinant+ minor
        
    return determinant








def circular_shift_matrix(matrix, shift_rows, shift_cols):
    """
    Performs a circular shift on a matrix.

    Args:
        matrix: The input matrix (NumPy array).
        shift_rows: The number of positions to shift rows (positive for down, negative for up).
        shift_cols: The number of positions to shift columns (positive for right, negative for left).

    Returns:
        A new matrix with the elements circularly shifted.
    """
    shifted_matrix = np.roll(matrix, shift=(shift_rows, -shift_cols), axis=(0, 1))
    return shifted_matrix





 #parallel speeds up computation only over very large matrices
# M is a mxn matrix binary matrix 
# all elements in M should be uint8 
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




def lift(matrix,G):
    
    num_rows = matrix.shape[0]
    
    num_cols=matrix.shape[1]
    
    liftSize=G
    
    output=np.zeros((liftSize*num_rows,liftSize*num_cols))
    
    
    for i in range(num_rows):
        
        for j in range(num_cols):
            
            circ=np.identity(liftSize)
            
            temp=np.zeros(liftSize)
            
            powers=matrix[i,j].nonzero_degrees
            
            
            coe=matrix[i,j].nonzero_coeffs
            
            for k in range(powers.shape[0]):
                
                temp=(temp+int(coe[k])*circular_shift_matrix(circ.copy(),0,powers[k]))%2
                
            output[i*liftSize:i*liftSize+liftSize,j*liftSize:j*liftSize+liftSize]=temp
            
    
    return output



def singlelift(poly,G):
    liftSize=G
    
    
    circ=np.identity(liftSize)
    
    temp=np.zeros((liftSize,liftSize))
    
    powers=poly.nonzero_degrees
    
    
    coe=poly.nonzero_coeffs
            
    for k in range(powers.shape[0]):
                
        temp=(temp+int(coe[k])*circular_shift_matrix(circ.copy(),0,powers[k]))%2
            
    
    return temp

    
    
def singleconjugate(poly,G):
    
    liftSize=G
    
    poly=singlelift(poly,G).transpose()
  
    
    output=[]
    
    basis=np.zeros((liftSize**2,liftSize))
    
    for i in range(np.shape(basis)[1]):
        
        eye=np.identity(liftSize)
        
        circ=circular_shift_matrix(eye.copy(), 0, i)
        
        basis[:,i]=circ.reshape(-1)
        
    A=basis.copy()
    
    b=poly.reshape(-1,1)
    
   
    A= np.append(A, b, axis=1)
    
    A=A.astype(int)
    
    x=gf2elim(A)[:liftSize, -1]
    
    x=x[::-1]
    
    return galois.Poly(x.transpose(),field=galois.GF(2))

    
    

def conjugate(matrix,G):
    
    num_rows = matrix.shape[0]
    
    num_cols=matrix.shape[1]
    
    liftSize=G
    
    matrix=lift(matrix,G).transpose()
    
    output=[]
    
    basis=np.zeros((liftSize**2,liftSize))
    
    for i in range(np.shape(basis)[1]):
        
        eye=np.identity(liftSize)
        
        circ=circular_shift_matrix(eye.copy(), 0, i)
        
        basis[:,i]=circ.reshape(-1)
 
    for i in range(num_cols):
        
        for j in range(num_rows):
            
            b=matrix[i*liftSize:i*liftSize+liftSize,j*liftSize:j*liftSize+liftSize].reshape(-1,1)
            
            A=basis.copy()
            
            A= np.append(A, b, axis=1)
            
            A=A.astype(int)
            
            x=gf2elim(A)[:liftSize, -1]
            
            x=x[::-1]
            
            output.append(galois.Poly(x.transpose(),field=galois.GF(2)))
            


    output=np.array(output)
    
    output=output.reshape(num_cols,num_rows)
    
    return output
    
    
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


