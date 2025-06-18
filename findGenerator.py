import galois
import numpy as np
import itertools
from itertools import permutations


def findGenerator(H,G):

  GF = galois.GF(2)
    
  num_rows = H.shape[0]
  
  num_cols=H.shape[1]

  Gen_rows=num_cols-num_rows
    
  assert num_rows>1, "!!!!Generator Matrix must Have at Least 2 Rows!!!!"
  assert Gen_rows>0, "!!!!Parity Check Matrix Full (Base) Rank!!!!"

  seed=list(range(num_cols))

  Generator=np.empty((Gen_rows,num_cols),dtype=object)
  
  
  all_combos = np.array(list(itertools.combinations(seed,num_rows+1)))
 
  for i in  range(Gen_rows):

    for j in range(num_cols):

        if j not in all_combos[i]:
        
            Generator[i,j]=galois.Poly([0],field=GF)
           
            
            break
        
        partialH=np.delete(H.copy(), j, axis=1)
        
        
        ciT=getMatrixDeterminant(partialH)
        
        ci=singleconjugate(ciT,G)
        
        
        Generator[i,j]=ci
       
 
  binaryGen=lift(Generator,G)

  extra_rows=G*Gen_rows-rankMat(binaryGen)
    
  nonzero_piece= [1] * G
    

  p1=galois.Poly(nonzero_piece,field=GF)
    
  p0=galois.Poly([0],field=GF)
    
  extrarow=[]
    
  for i in range(num_cols):
      
      if i==0 or i==1:
          
          extrarow.append(p1)
          
      else:
          
          extrarow.append(p0)
  
  extrarow=np.array(extrarow)
    
  perms = np.array(list(set(permutations(extrarow))))
  assert extra_rows<perms.shape[0], "Base Matrix has too many linearly dependent rows over Polynomial Ring." 
          
  i=0 
    
  while extra_rows>0:
      
      Generator=np.insert(Generator,Generator.shape[0] ,perms[i,:],axis= 0)
      
      binaryGen=lift(Generator,G)
      
      extra_rows=G*Gen_rows-rankMat(binaryGen)
      i+=1
  


  return Generator
        
          

  
      
      
    
    
    
    
