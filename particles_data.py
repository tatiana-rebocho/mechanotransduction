import random

for i in range(1,17):

  n = 1200 #number of particles
    
  #data file name, where i is the bed index
  file_name = str(i)+ ".bed" + ".data" 

  #open file file_name to write it
  f = open(file_name, "w+")

  tp = 1 #number of particle types

  #printing of information needed in data file
  f.write("# lammps molecular data \n \n")

  #definition of number and types of particles and bonds
  f.write(str(n) + " atoms \n \n")
  f.write(str(0) + " bonds \n \n")

  f.write(str(tp) + " atom types \n \n")
  f.write(str(0) + " bond types \n \n")
    
  #definition of the mean diameter
  diameter = 1.00
    
  #definition of the dimensions of the box
  dim = 20 * diameter

  f.write("0 "+str(dim)+" xlo xhi \n") #dimensions in x
  f.write("0 "+str(dim)+" ylo yhi \n") #dimensions in y
  f.write("0 "+str(dim)+" zlo zhi \n \n \n") #dimensions in z

  f.write("Atoms \n \n")
    
  #initializes the molecule counting (each particle will be part of a molecule)
  mol = 1

  for k in range (1, n+1): #cycle for particle ID
    x = round(random.uniform(0,dim), 4) #position in x
    y = round(random.uniform(0,dim), 4) #position in y
    z = round(random.uniform(2/20*dim,10/20*dim), 4) #position in z
    diam = random.gauss(diameter, 0.05) #particle diameter
    den = 1.000 #particle density
        
    #prints a line per particle 
    print (k, 1, x, y, z, mol, diam, den, file=f)
    #particleID particleType x y z moleculeID diameter density

    #increases molecule ID
    mol = mol + 1
        
  f.close()