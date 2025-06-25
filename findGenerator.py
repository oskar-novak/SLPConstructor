
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
           
            
            continue
        temp=all_combos[i].copy()
        temp=temp[temp != j]
        
        
        partialH=H.copy()
        partialH=partialH[:,temp]
        
        
        
        
        ciT=getMatrixDeterminant(partialH)
        
      
        
        ci=singleconjugate(ciT,G)
        
        
        Generator[i,j]=ci
       
  
  binaryGen=lift(Generator,G)
  liftH=lift(H,G)

  extra_rows=(liftH.shape[1]-rankMat(liftH))-rankMat(binaryGen)
    
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
    
  
  good_rows=[]
  
  for i in range(perms.shape[0]):
      
      
      checkArray=lift(H@conjugate(np.asmatrix(perms[i,:]),G),G)
      if np.all(checkArray==0):
         good_rows.append(i)

  if extra_rows<=len(good_rows):
      
      Generator=np.insert(Generator,Generator.shape[0] ,perms[good_rows[:extra_rows],:],axis= 0)
      
  else:
      
      Generator=np.insert(Generator,Generator.shape[0] ,perms[good_rows,:],axis= 0)
      
    
    
  binaryGen=lift(Generator,G)
    
  if (liftH.shape[1]-rankMat(liftH))-rankMat(binaryGen)>0:
      print("!!!Found Sub Code Generator Matrix with Rank Deficiency {}!!!".format((liftH.shape[1]-rankMat(liftH))-rankMat(binaryGen)))
          
  
  return Generator
        
          

  
      
      
    
    
    
    
