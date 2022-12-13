import random
import math

for i in range(1,11):

  name = str(i) + ".cell.data" #data file name, where i is the cell index
       
  f = open(name, "w+") #open file file_name to write it

  n_atoms = 1 + 6 + 12 + 24 #number of elements per cell
  
  n_bonds = 12 + 24 + 48 #number of bonds per cell
    
  type_atoms = 4 #number of particles types (bed+1 per layer)
    
  bond_types = 2 + 2 + 2 #number of bond types
    
  #printing of information needed in data file
  f.write("# lammps molecular data \n \n")

  #definition of number and types of particles and bonds
  f.write(str(n_atoms) + " atoms \n")
  f.write(str(n_bonds) + " bonds \n \n")

  f.write(str(type_atoms) + " atom types \n")
  f.write(str(bond_types) + " bond types \n \n")

  #initializes the coordinates section
  f.write("Atoms \n \n")

  #type 1

  k = 1 #initializes the particle counting

  mol = 1 #defines the molecule counting (1 per cell)
    
  diam_mp = 1.00 #definition of the mean diameter for microparticles
  dim = 20 * diam_mp #definition of the dimensions of the box
    
  #definition of the limits where the cell can exist
    
  x_inflim = 6.5
  x_suplim = dim - 6.5
  y_inflim = 6.5
  y_suplim = dim - 6.5
    
  #randomly chooses the center of the cell (x_o,y_o)
  x_o = random.uniform(x_inflim,x_suplim)
  y_o = random.uniform(y_inflim,y_suplim)
    
  z = diam_mp * 5 #defines the at which height the cell will start z
    
  #definition of the diameter of each subcellular elements its radius and its density
  diam = 1.5
  radius = diam / 2
  den = 1.000
    
  #prints a line per particle 
  print (k, 2, x_o, y_o, z, mol, diam, den, file=f) 
  #atom-ID atom-type x y z moleculeID diameter density 
    
  #type 2
   
  n_layer = 6 #number of entities per layer
  #initializes the list at which angle the subcellular element will be
  teta_list = []
  #defines the first angle (in rad)
  teta = 0
    
  #for each element there will be an angle teta
  for n in range(n_layer):
    teta_list.append(teta)
    teta = teta + 2 * math.pi/ n_layer
    
  #for each angle we define the element position x and y
  for teta in teta_list: #particle ID
    x = x_o + radius * math.sin(teta)
    y = y_o + radius * math.cos(teta)
        
    #increases atom ID
    k = k + 1
        
    #prints a line per particle 
    print (k, 2, x, y, z, mol, diam, den, file=f)
    #atom-ID atom-type x y z moleculeID diameter density 
        
  #type 3

  n_layer = 12 #number of cells per layer
  #initializes the list at which angle the subcellular element will be
  teta_list = []
  #defines the first angle (in rad)
  teta = - math.pi / n_layer
    
  #for each element there will be an angle teta
  for n in range(n_layer):
    teta_list.append(teta)
    teta = teta + 2 * math.pi/ n_layer
    
  #for each angle we define the element position x and y
  for teta in teta_list: #particle ID
    x = x_o + 2*radius * math.sin(teta)
    y = y_o + 2*radius * math.cos(teta)
        
    #increases atom ID
    k = k + 1
        
    #prints a line per particle 
    print (k, 3, x, y, z, mol, diam, den, file=f) 
    #atom-ID atom-type x y z moleculeID diameter density 
        
  #type 4

  n_layer = 24 #number of cells per layer
  #initializes the list at which angle the subcellular element will be
  teta_list = []
  #defines the first angle (in rad)
  teta = - 3 * math.pi / n_layer
    
  #for each element there will be an angle teta
  for n in range(n_layer):
    teta_list.append(teta)
    teta = teta + 2 * math.pi/ n_layer
  
  #for each angle we define the element position x and y
  for teta in teta_list: #particle ID
    x = x_o + 3*radius * math.sin(teta)
    y = y_o + 3*radius * math.cos(teta)
        
    #increases atom ID
    k = k + 1
        
    #prints a line per particle 
    print (k, 4, x, y, z, mol, diam, den, file=f) 
    #atom-ID atom-type x y z moleculeID diameter density 
        
  #initializes the bonds section
        
  f.write("\n \n")
    
  f.write("Bonds \n \n")

  #define the atomID per layer
  type_1 = [1]
  type_2 = [*range(2,7+1)]
  type_3 = [*range(8,19+1)]
  type_4 = [*range(20,43+1)]

  #initializes the IDbond and bondType counting
  ID_bond = 1
  bond_type = 1
    
  #bonds for type1 with type2

  #for each element in type2 there will be a bond with the element in type1
  for x in type_2: #prints a line per bond
    print(ID_bond, bond_type, type_1[0], x, file=f)
    #bondID bondType elementID1 elementID2
    ID_bond = ID_bond + 1 #increases bond ID

  #bonds for type2 with type2

  bond_type = bond_type + 1 #increases bond type
    
  #for each element in type2 there will be a bond with the two elements next to it in the same layer
  for n in range(0,len(type_2)-1):
    #for the first element in the list it will be bonded with the last element in it and with the next one
    if n == 0:
      print(ID_bond, bond_type, type_2[n], type_2[len(type_2)-1], file=f)
      #bondID bondType elementID1 elementID2
      ID_bond = ID_bond + 1 #increases bond ID
      print(ID_bond, bond_type, type_2[n], type_2[n+1], file=f)
      #bondID bondType elementID1 elementID2
      ID_bond = ID_bond + 1 #increases bond ID
    else:
      print(ID_bond, bond_type, type_2[n], type_2[n+1], file=f)
      #bondID bondType elementID1 elementID2
      ID_bond = ID_bond + 1 #increases bond ID

  #bonds for type2 with type3

  bond_type = bond_type + 1 #increases bond type

  k = type_3[0] #k is the first element ID for list type3

  #while k is less or equal to the last element ID of type3
  while k <= type_3[len(type_3)-1]:
    for x in type_2: #for each element in type2
      #two elements of type3 will be bonded to an element of type2
      #j counts the two bonds
      j = 1
      while j <=2:
        print(ID_bond, bond_type, x, k, file=f)
        #bondID bondType elementID1 elementID2
        ID_bond = ID_bond + 1 #increases bond ID
        #increases k and j
        k = k + 1
        j = j + 1

  #bonds for type3 with type3

  bond_type = bond_type + 1 #increases bond type

  #for each element in type3 there will be a bond with the two elements next to it in the same layer
  for n in range(0,len(type_3)-1):
    #for the first element in the list it will be bonded with the last element in it and with the next one
    if n == 0:
      print(ID_bond, bond_type, type_3[n], type_3[len(type_3)-1], file=f)
      #bondID bondType elementID1 elementID2
      ID_bond = ID_bond + 1 #increases bond ID
      print(ID_bond, bond_type, type_3[n], type_3[n+1], file=f)
      #bondID bondType elementID1 elementID2
      ID_bond = ID_bond + 1 #increases bond ID
    else:
      print(ID_bond, bond_type, type_3[n], type_3[n+1], file=f)
      #bondID bondType elementID1 elementID2
      ID_bond = ID_bond + 1 #increases bond ID

  #bonds for type3 with type4

  bond_type = bond_type + 1 #increases bond type

  k = type_4[0] #k is the first element ID for list type4
    
  #while k is less or equal to the last element ID of type4
  while k <= type_4[len(type_4)-1]:
    for x in type_3: #for each element in type3
      #two elements of type3 will be bonded to an element of type2
      #j counts the two bonds
      j = 1
      while j <=2:
        print(ID_bond, bond_type, x, k, file=f)
        #bondID bondType elementID1 elementID2
        ID_bond = ID_bond + 1 #increases bond ID
        #increases k and j
        k = k + 1
        j = j + 1

  #bonds for type4 with type4

  bond_type = bond_type + 1 #increases bond type

  #for each element in type4 there will be a bond with the two elements next to it in the same layer
  for n in range(0,len(type_4)-1):
    #for the first element in the list it will be bonded with the last element in it and with the next one
    if n == 0:
      print(ID_bond, bond_type, type_4[n], type_4[len(type_4)-1], file=f)
      #bondID bondType elementID1 elementID2
      ID_bond = ID_bond + 1 #increases bond ID
      print(ID_bond, bond_type, type_4[n], type_4[n+1], file=f)
      #bondID bondType elementID1 elementID2
      ID_bond = ID_bond + 1 #increases bond ID
    else:
      print(ID_bond, bond_type, type_4[n], type_4[n+1], file=f)
      #bondID bondType elementID1 elementID2
      ID_bond = ID_bond + 1 #increases bond ID
    
  f.close()
