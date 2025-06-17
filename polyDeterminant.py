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
