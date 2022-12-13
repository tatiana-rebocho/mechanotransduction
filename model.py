import shutil
import math
import os.path

from lammps import lammps

thr_freeze = 10**(-4)

size_dim = SIZE
bed_ID = BEDID
runidx = [ID_FILE]

for run in runidx:

  lmp = lammps()

  lmp.file("variables.run")

  lmp.command("variable bedID equal "+str(bed_ID))
  lmp.command("variable runindx equal "+str(run))
  lmp.command("variable size equal "+str(size_dim))
    
  lmp.file("initialization.run")
  lmp.file("bed-definition.run")
  lmp.file("cell-inclusion.run")
   
  #create patches
    
  cell = range(1201,1243+1)
  active_layer = range(1201,1207+1)

  nlocal = lmp.extract_global("nlocal")
  ids = lmp.extract_atom("id")
  ids_list = []
    
  for index in range(nlocal):
    ids_list.append(ids[index]) 
        
  n_atoms = max(ids_list)    
    
  pair_cell_patch = []
  patches_IDs = []
  mp_patch =[]
    
  for i in range(0,600):
    nlocal = lmp.extract_global("nlocal")
    x = lmp.extract_atom("x")
    ids = lmp.extract_atom("id")
    rad = lmp.extract_atom("radius")
    ids_list = []
        
    for index in range(nlocal):
      ids_list.append(ids[index])
            
    for IDmp in range(1,1201):
      i_mp = ids_list.index(IDmp)
            
        for ID in range(1201,1244):
          i_cell = ids_list.index(ID)   
          Xmp = x[i_mp][0]
          X = x[i_cell][0]
          Ymp = x[i_mp][1]
          Y = x[i_cell][1]
          Zmp = x[i_mp][2]
          Z= x[i_cell][2]    
          DX = (Xmp-X)**2.0
          DY = (Ymp-Y)**2.0
          DZ = (Zmp-Z)**2.0 
          dist = (DX+DY+DZ)**(0.5)
          rad_sum = rad[i_cell] + rad[i_mp]
          pair = [IDmp,ID]
                
          if dist <= rad_sum and pair not in mp_patch:
            mp_patch.append(pair)
            frac = rad[i_cell]/ dist     
            vec_x = x[i_mp][0]-x[i_cell][0]
            vec_y = x[i_mp][1]-x[i_cell][1]
            vec_z = x[i_mp][2]-x[i_cell][2]    
            x_p = x[i_cell][0] + frac * vec_x
            y_p = x[i_cell][1] + frac * vec_y
            z_p = x[i_cell][2] + frac * vec_z   
            lmp.command("create_atoms 8 single " +str(x_p)+ " "+str(y_p)+ " "+str(z_p))
            n_atoms = n_atoms + 1
            dist_patch = ((x[i_cell][0]-x_p)**2.0+(x[i_cell][1]-y_p)**2.0+(x[i_cell][2]-z_p)**2.0)**0.5
            lmp.command("set atom "+str(n_atoms)+ " mol " + str(IDmp))     
            pair_cell_patch.append([ID,n_atoms,dist_patch])
            patches_IDs.append(n_atoms)
                    
            if ID in active_layer:
              lmp.command("create_bonds single/bond 7 " + str(ID) + " " + str(n_atoms))
                        
            else:
              lmp.command("create_bonds single/bond 8 " + str(ID) + " " + str(n_atoms))
    
    lmp.command("set type 8 diameter 0.1")
    lmp.command("set type 8 mass 0.000001")   
    lmp.command("group MPpatches type 8")
    lmp.command("group FreezeNow union FreezeNow MPpatches")
    lmp.command("group MPpatches include molecule")
    lmp.command("group langevin subtract all MPpatches")
    lmp.command("fix integrator langevin nve")
    lmp.command("fix 1 langevin langevin ${temperature} ${temperature} ${damp} ${seed}")
    lmp.command("fix patches MPpatches rigid/small molecule langevin ${temperature} ${temperature} ${damp} ${seed}")
    
    lmp.command("run 100")

  #correct patches

  lmp.command("group del type 8")
  lmp.command("delete_atoms group del compress no bond yes")
  lmp.command("group del delete")
    
  nlocal = lmp.extract_global("nlocal")
  x = lmp.extract_atom("x")
  ids = lmp.extract_atom("id")
  rad = lmp.extract_atom("radius")
  ids_list = []
    
  for index in range(nlocal):
    ids_list.append(ids[index])
    
  n_atoms = max(ids_list)
  
  cell_patch = []
  patches_IDs = []
  distance_patchs = []
    
  for IDmp in range(1,1201):
    i_mp = ids_list.index(IDmp)
        
    for ID in range(1201,1244):
        
      if ID in cell:
        i_cell = ids_list.index(ID)  
        Xmp = x[i_mp][0]
        X = x[i_cell][0]
        Ymp = x[i_mp][1]
        Y = x[i_cell][1]
        Zmp = x[i_mp][2]
        Z= x[i_cell][2]
        DX = (Xmp-X)**2.0
        DY = (Ymp-Y)**2.0
        DZ = (Zmp-Z)**2.0
        dist = (DX+DY+DZ)**(0.5)
        rad_sum = rad[i_cell] + rad[i_mp]
                
        if dist <= rad_sum:
          frac = rad[i_cell]/ dist     
          vec_x = x[i_mp][0]-x[i_cell][0]
          vec_y = x[i_mp][1]-x[i_cell][1]
          vec_z = x[i_mp][2]-x[i_cell][2]  
          x_p = x[i_cell][0] + frac * vec_x
          y_p = x[i_cell][1] + frac * vec_y
          z_p = x[i_cell][2] + frac * vec_z
          lmp.command("create_atoms 8 single " +str(x_p)+ " "+str(y_p)+ " "+str(z_p))
          n_atoms = n_atoms + 1
          dist_patch = ((x[i_cell][0]-x_p)**2.0+(x[i_cell][1]-y_p)**2.0+(x[i_cell][2]-z_p)**2.0)**0.5
          lmp.command("set atom "+str(n_atoms)+ " mol " + str(IDmp))       
          cell_patch.append(ID)
          patches_IDs.append(n_atoms)
          distance_patchs.append(dist_patch)
                    
          if ID in active_layer:
            lmp.command("create_bonds single/bond 7 " + str(ID) + " " + str(n_atoms))
                        
          else:
            lmp.command("create_bonds single/bond 8 " + str(ID) + " " + str(n_atoms))
    
  lmp.command("set type 8 diameter 0.1")
  lmp.command("set type 8 mass 0.000001")
  lmp.command("unfix integrator")
  lmp.command("unfix 1")
  lmp.command("unfix patches")
  lmp.command("group MPpatches delete")
  lmp.command("group langevin delete")
  lmp.command("group MPpatches type 8")
  lmp.command("group FreezeNow union FreezeNow MPpatches")
  lmp.command("group MPpatches include molecule")
  lmp.command("group langevin subtract all MPpatches")
  lmp.command("fix integrator langevin nve")
  lmp.command("fix 1 langevin langevin ${temperature} ${temperature} ${damp} ${seed}")
  lmp.command("fix patches MPpatches rigid/small molecule langevin ${temperature} ${temperature} ${damp} ${seed}")
    
  lmp.command("run 100")
    
  #contraction 1
    
  #tree bonds
    
  r1 = lmp.extract_variable("r1",group=None, vartype=0)
  r2 = lmp.extract_variable("r2",group=None, vartype=0)
  r3 = lmp.extract_variable("r3",group=None, vartype=0)
  r4 = lmp.extract_variable("r4",group=None, vartype=0)
  r5 = lmp.extract_variable("r5",group=None, vartype=0)
  r6 = lmp.extract_variable("r6",group=None, vartype=0)

  centerID = 1201
  layer1_ID = range(1202,1207+1)
  layer2_ID = range(1208,1219+1)
  layer3_ID = range(1220,1243+1)
    
  layer1 = range(1202,1207+1)
  layer2 = range(1208,1219+1)
  layer3 = range(1220,1243+1)

  k = layer2[0]
  q = layer3[0]

  tree_bonds = {i:[] for i in layer1_ID}

  for ID in layer1_ID:
    tree_bonds[ID].append(ID)
    j = 1
    
    while j<=2:
      tree_bonds[ID].append(k)
      k = k + 1
      j = j + 1
    l = 1
    
    while l<=4:
      tree_bonds[ID].append(q)
      q = q + 1
      l = l + 1
        
  #calculate distances 
  distancesBegin = []
  IDBegin = []
    
  for ID in layer1_ID:
    i_c = ids_list.index(centerID)
    i_n = ids_list.index(ID)
    centerX = x[i_c][0]
    X = x[i_n][0]
    centerY = x[i_c][1]
    Y = x[i_n][1]
    centerZ = x[i_c][2]
    Z= x[i_n][2]
    DX = (centerX-X)**2.0
    DY = (centerY-Y)**2.0
    DZ = (centerZ-Z)**2.0
    dist = (DX+DY+DZ)**(0.5)
    distancesBegin.append(dist)  
    IDBegin.append(ID)
    
  desactive = []
      
  for i in range(0,10000):
    dt = lmp.get_thermo("dt")
    dr1 = r1-r1*0.1*i*dt
    dr2 = r2-r2*0.1*i*dt
    lmp.command("bond_style harmonic")
    lmp.command("bond_coeff 1 ${bondK} "+str(dr1))
    lmp.command("bond_coeff 2 ${bondK} "+str(dr2))
    lmp.command("bond_coeff 3 ${bondKdr} ${r3}")
    lmp.command("bond_coeff 4 ${bondKdr} ${r4}")
    lmp.command("bond_coeff 5 ${bondKdr} ${r5}")
    lmp.command("bond_coeff 6 ${bondKdr} ${r6}")
    lmp.command("bond_coeff 7 ${bondmc} ${r7}")
    lmp.command("bond_coeff 8 ${bondinit} ${r7}")
    nlocal = lmp.extract_global("nlocal")
    x = lmp.extract_atom("x")
    ids = lmp.extract_atom("id")
    ids_list = []
        
    for index in range(nlocal):
      ids_list.append(ids[index])
        
    for ID in range(1201,1244):
        
      if ID in active_layer:
            
        if ID in cell_patch:
          idx_cell = cell_patch.index(ID)
          i_cell = ids_list.index(ID)
          i_p = ids_list.index(patches_IDs[idx_cell])
          d0 = distance_patchs[idx_cell]
          cellX = x[i_cell][0]
          X = x[i_p][0]
          cellY = x[i_cell][1]
          Y = x[i_p][1]
          cellZ = x[i_cell][2]
          Z= x[i_p][2]      
          DX = (cellX-X)**2.0
          DY = (cellY-Y)**2.0
          DZ = (cellZ-Z)**2.0
          dist = (DX+DY+DZ)**(0.5)     
          diff = abs(d0-dist)
                    
          if diff >= thr_freeze:
            lmp.command("set atom " + str(ID)+" type 7")
            lmp.command("group zeroforce type 7")
            lmp.command("group langevin subtract all zeroforce MPpatches")
                        
            if ID in active_layer:
              active_layer.remove(ID)
                        
        if ID in IDBegin:
                
          if ID in active_layer:
            i_begin = IDBegin.index(ID)
            i_c = ids_list.index(centerID)
            i_n = ids_list.index(ID)          
            centerX = x[i_c][0]
            X = x[i_n][0]
            centerY = x[i_c][1]
            Y = x[i_n][1]
            centerZ = x[i_c][2]
            Z= x[i_n][2]       
            DX = (centerX-X)**2.0
            DY = (centerY-Y)**2.0
            DZ = (centerZ-Z)**2.0      
            dist = (DX+DY+DZ)**(0.5)      
            thr = distancesBegin[i_begin] * 0.96
                        
            if dist < thr:
              lmp.command("set atom "+str(ID)+" type 6")
              active_layer.remove(ID)
              layer1_ID.remove(ID)
              del_atoms = tree_bonds[ID]
                            
              if ID not in desactive:
                desactive.append(ID)
                                
              for new_ID in del_atoms:
                lmp.command("set atom "+str(new_ID)+" type 6")
                                
                if new_ID not in desactive:
                  desactive.append(new_ID)
                 
    lmp.command("group del type 6")
    lmp.command("delete_atoms group del compress no bond yes")
    lmp.command("group del delete")
                 
    lmp.command("run 5")

  lmp.command("variable dr1 equal "+str(dr1))
  lmp.command("variable dr2 equal "+str(dr2))    

  #next layer

  lmp.command("group del type 8")
  lmp.command("delete_atoms group del compress no bond yes")
  lmp.command("group del delete")
    
  lmp.file("next-layer-1.run")
    
  active_layer = range(1208,1219+1)

  for ID in desactive:
    if ID in active_layer:
      active_layer.remove(ID)
            
  cell = range(1201,1243+1)

  for ID in desactive:
    if ID in cell:
      cell.remove(ID)

  #create patches
    
  nlocal = lmp.extract_global("nlocal")
  ids = lmp.extract_atom("id")
  ids_list = []
    
  for index in range(nlocal):
    ids_list.append(ids[index])
    
  n_atoms = max(ids_list)
    
  pair_cell_patch = []
  patches_IDs = []
  mp_patch =[]
    
  for i in range(0,300):
    nlocal = lmp.extract_global("nlocal")
    x = lmp.extract_atom("x")
    ids = lmp.extract_atom("id")
    rad = lmp.extract_atom("radius")
    ids_list = []
        
    for index in range(nlocal):
      ids_list.append(ids[index])
    
    for IDmp in range(1,1201):
      i_mp = ids_list.index(IDmp)
            
      for ID in range(1201,1244):
            
        if ID in cell:
          i_cell = ids_list.index(ID)        
          Xmp = x[i_mp][0]
          X = x[i_cell][0]
          Ymp = x[i_mp][1]
          Y = x[i_cell][1]
          Zmp = x[i_mp][2]
          Z= x[i_cell][2]      
          DX = (Xmp-X)**2.0
          DY = (Ymp-Y)**2.0
          DZ = (Zmp-Z)**2.0   
          dist = (DX+DY+DZ)**(0.5)
          rad_sum = rad[i_cell] + rad[i_mp]
          pair = [IDmp,ID]
                    
          if dist <= rad_sum and pair not in mp_patch:
            mp_patch.append(pair)
            frac = rad[i_cell]/ dist      
            vec_x = x[i_mp][0]-x[i_cell][0]
            vec_y = x[i_mp][1]-x[i_cell][1]
            vec_z = x[i_mp][2]-x[i_cell][2]      
            x_p = x[i_cell][0] + frac * vec_x
            y_p = x[i_cell][1] + frac * vec_y
            z_p = x[i_cell][2] + frac * vec_z     
            lmp.command("create_atoms 8 single " +str(x_p)+ " "+str(y_p)+ " "+str(z_p))
            n_atoms = n_atoms + 1
            dist_patch = ((x[i_cell][0]-x_p)**2.0+(x[i_cell][1]-y_p)**2.0+(x[i_cell][2]-z_p)**2.0)**0.5
            lmp.command("set atom "+str(n_atoms)+ " mol " + str(IDmp))      
            pair_cell_patch.append([ID,n_atoms,dist_patch])
            patches_IDs.append(n_atoms)
                        
            if ID in active_layer or ID in layer1 or ID==centerID:
              lmp.command("create_bonds single/bond 7 " + str(ID) + " " + str(n_atoms))
                            
            else:
              lmp.command("create_bonds single/bond 8 " + str(ID) + " " + str(n_atoms))
                        
    lmp.command("set type 8 diameter 0.1")
    lmp.command("set type 8 mass 0.000001")
    lmp.command("unfix integrator")
    lmp.command("unfix 1")
    lmp.command("unfix patches")
    lmp.command("group MPpatches delete")
    lmp.command("group langevin delete")
    lmp.command("group MPpatches type 8")
    lmp.command("group FreezeNow union FreezeNow MPpatches")
    lmp.command("group MPpatches include molecule")
    lmp.command("group langevin subtract all MPpatches zeroforce")
    lmp.command("fix integrator langevin nve")
    lmp.command("fix 1 langevin langevin ${temperature} ${temperature} ${damp} ${seed}")
    lmp.command("fix patches MPpatches rigid/small molecule langevin ${temperature} ${temperature} ${damp} ${seed}")
        
    lmp.command("run 100")
        
  #correct patches

  lmp.command("group del type 8")
  lmp.command("delete_atoms group del compress no bond yes")
  lmp.command("group del delete")
  nlocal = lmp.extract_global("nlocal")
  x = lmp.extract_atom("x")
  ids = lmp.extract_atom("id")
  rad = lmp.extract_atom("radius")
  ids_list = []
    
  for index in range(nlocal):
    ids_list.append(ids[index])
    
  n_atoms = max(ids_list)
    
  cell_patch = []
  patches_IDs = []
  distance_patchs = []
    
  for IDmp in range(1,1201):
    i_mp = ids_list.index(IDmp)
        
    for ID in range(1201,1244):
        
      if ID in cell:
        i_cell = ids_list.index(ID)   
        Xmp = x[i_mp][0]
        X = x[i_cell][0]
        Ymp = x[i_mp][1]
        Y = x[i_cell][1]
        Zmp = x[i_mp][2]
        Z= x[i_cell][2]
        DX = (Xmp-X)**2.0
        DY = (Ymp-Y)**2.0
        DZ = (Zmp-Z)**2.0     
        dist = (DX+DY+DZ)**(0.5)
        rad_sum = rad[i_cell] + rad[i_mp]
                
        if dist <= rad_sum:
          frac = rad[i_cell]/ dist
          vec_x = x[i_mp][0]-x[i_cell][0]
          vec_y = x[i_mp][1]-x[i_cell][1]
          vec_z = x[i_mp][2]-x[i_cell][2]      
          x_p = x[i_cell][0] + frac * vec_x
          y_p = x[i_cell][1] + frac * vec_y
          z_p = x[i_cell][2] + frac * vec_z     
          lmp.command("create_atoms 8 single " +str(x_p)+ " "+str(y_p)+ " "+str(z_p))
          n_atoms = n_atoms + 1
          dist_patch = ((x[i_cell][0]-x_p)**2.0+(x[i_cell][1]-y_p)**2.0+(x[i_cell][2]-z_p)**2.0)**0.5
          lmp.command("set atom "+str(n_atoms)+ " mol " + str(IDmp))         
          cell_patch.append(ID)
          patches_IDs.append(n_atoms)
          distance_patchs.append(dist_patch)
                    
          if ID in active_layer or ID in layer1 or ID==centerID:
            lmp.command("create_bonds single/bond 7 " + str(ID) + " " + str(n_atoms))
                        
          else:
            lmp.command("create_bonds single/bond 8 " + str(ID) + " " + str(n_atoms))
                    
  lmp.command("set type 8 diameter 0.1")
  lmp.command("set type 8 mass 0.000001")
  lmp.command("unfix integrator")
  lmp.command("unfix 1")
  lmp.command("unfix patches")
  lmp.command("group MPpatches delete")
  lmp.command("group langevin delete")
  lmp.command("group MPpatches type 8")
  lmp.command("group FreezeNow union FreezeNow MPpatches")
  lmp.command("group MPpatches include molecule")
  lmp.command("group langevin subtract all MPpatches zeroforce")
  lmp.command("fix integrator langevin nve")
  lmp.command("fix 1 langevin langevin ${temperature} ${temperature} ${damp} ${seed}")
  lmp.command("fix patches MPpatches rigid/small molecule langevin ${temperature} ${temperature} ${damp} ${seed}")
    
  lmp.command("run 100")
    
  #contraction 2
    
  #tree bonds
    
  k = layer2[0]
  q = layer3[0]

  tree_bonds = {i:[] for i in layer2_ID}
  
  for ID in layer2_ID:
    tree_bonds[ID].append(ID)
    l = 1
        
    while l<=2:
      tree_bonds[ID].append(q)
      q = q + 1
      l = l + 1

  for ID in desactive:
    if ID in layer2_ID:
      layer2_ID.remove(ID)
    
  #calculate distances    
  distancesBegin = []
  IDBegin = []
    
  if len(layer2_ID) != 0:
    
    for ID in layer2_ID:
      i_c = ids_list.index(centerID)
      i_n = ids_list.index(ID)
      centerX = x[i_c][0]
      X = x[i_n][0]
      centerY = x[i_c][1]
      Y = x[i_n][1]
      centerZ = x[i_c][2]
      Z= x[i_n][2]   
      DX = (centerX-X)**2.0
      DY = (centerY-Y)**2.0
      DZ = (centerZ-Z)**2.0  
      dist = (DX+DY+DZ)**(0.5)
      distancesBegin.append(dist)     
      IDBegin.append(ID)

  for i in range(0,10000):
    dt = lmp.get_thermo("dt")
    dr3 = r3-r3*0.1*i*dt
    dr4 = r4-r4*0.1*i*dt
    lmp.command("bond_style harmonic")
    lmp.command("bond_coeff 1 ${bondK} "+str(dr1))
    lmp.command("bond_coeff 2 ${bondK} "+str(dr2))
    lmp.command("bond_coeff 3 ${bondK} "+str(dr3))
    lmp.command("bond_coeff 4 ${bondK} "+str(dr4))
    lmp.command("bond_coeff 5 ${bondKdr} ${r5}")
    lmp.command("bond_coeff 6 ${bondKdr} ${r6}")
    lmp.command("bond_coeff 7 ${bondmc} ${r7}")
    lmp.command("bond_coeff 8 ${bondinit} ${r7}")
    nlocal = lmp.extract_global("nlocal")
    x = lmp.extract_atom("x")
    ids = lmp.extract_atom("id")
    ids_list = []
        
    for index in range(nlocal):
      ids_list.append(ids[index])
            
    for ID in range(1201,1244):
        
      if ID in active_layer:
            
        if ID in cell_patch:
          idx_cell = cell_patch.index(ID)
          i_cell = ids_list.index(ID)
          i_p = ids_list.index(patches_IDs[idx_cell])
          d0 = distance_patchs[idx_cell]
          cellX = x[i_cell][0]
          X = x[i_p][0]
          cellY = x[i_cell][1]
          Y = x[i_p][1]
          cellZ = x[i_cell][2]
          Z= x[i_p][2]          
          DX = (cellX-X)**2.0
          DY = (cellY-Y)**2.0
          DZ = (cellZ-Z)**2.0
          dist = (DX+DY+DZ)**(0.5)       
          diff = abs(d0-dist)
                    
          if diff >= thr_freeze:
            lmp.command("set atom " + str(ID)+" type 7")
            lmp.command("group zeroforce type 7")
            lmp.command("group langevin subtract all zeroforce MPpatches")
                        
            if ID in active_layer:
              active_layer.remove(ID)

        if ID in active_layer:
                
          if ID in IDBegin:
            i_begin = IDBegin.index(ID)
            i_c = ids_list.index(centerID)
            i_n = ids_list.index(ID)       
            centerX = x[i_c][0]
            X = x[i_n][0]
            centerY = x[i_c][1]
            Y = x[i_n][1]
            centerZ = x[i_c][2]
            Z= x[i_n][2]
            DX = (centerX-X)**2.0
            DY = (centerY-Y)**2.0
            DZ = (centerZ-Z)**2.0        
            dist = (DX+DY+DZ)**(0.5)           
            thr = distancesBegin[i_begin] * 0.96
                        
            if dist < thr:
              lmp.command("set atom "+str(ID)+" type 6")
              active_layer.remove(ID)
              layer2_ID.remove(ID)
              del_atoms = tree_bonds[ID]
                            
              if ID not in desactive:
                desactive.append(ID)
                                
              for new_ID in del_atoms:
                lmp.command("set atom "+str(new_ID)+" type 6")
                                
                if new_ID not in desactive:
                  desactive.append(new_ID)
                        
    lmp.command("group del type 6")
    lmp.command("delete_atoms group del compress no bond yes")
    lmp.command("group del delete")
        
    lmp.command("run 5")
    
  lmp.command("variable dr3 equal "+str(dr3))
  lmp.command("variable dr4 equal "+str(dr4)) 

  ids = lmp.extract_atom("id")
  nlocal = lmp.extract_global("nlocal")
  ids_list = []
    
  for index in range(nlocal):
    ids_list.append(ids[index]) 
        
  #next layer
    
  lmp.command("group del type 8")
  lmp.command("delete_atoms group del compress no bond yes")
  lmp.command("group del delete")

  lmp.file("next-layer-2.run")

  active_layer = range(1220,1243+1)
    
  for ID in desactive:
    if ID in active_layer:
      active_layer.remove(ID)
    
  cell = range(1201,1243+1)

  for ID in desactive:
    if ID in cell:
      cell.remove(ID)
    
  #create patches
    
  nlocal = lmp.extract_global("nlocal")
  ids = lmp.extract_atom("id")
  ids_list = []
    
  for index in range(nlocal):
    ids_list.append(ids[index])
   
  n_atoms = max(ids_list)
    
  pair_cell_patch = []
  patches_IDs = []
    
  mp_patch =[]
    
  for i in range(0,300):
    nlocal = lmp.extract_global("nlocal")
    x = lmp.extract_atom("x")
    ids = lmp.extract_atom("id")
    rad = lmp.extract_atom("radius")
    ids_list = []
        
    for index in range(nlocal):
      ids_list.append(ids[index])
    
    for IDmp in range(1,1201):
      i_mp = ids_list.index(IDmp)
            
      for ID in range(1201,1244):
            
        if ID in cell:
          i_cell = ids_list.index(ID)    
          Xmp = x[i_mp][0]
          X = x[i_cell][0]
          Ymp = x[i_mp][1]
          Y = x[i_cell][1]
          Zmp = x[i_mp][2]
          Z= x[i_cell][2]    
          DX = (Xmp-X)**2.0
          DY = (Ymp-Y)**2.0
          DZ = (Zmp-Z)**2.0      
          dist = (DX+DY+DZ)**(0.5)
          rad_sum = rad[i_cell] + rad[i_mp]
          pair = [IDmp,ID]
                    
          if dist <= rad_sum and pair not in mp_patch:
            mp_patch.append(pair)
            frac = rad[i_cell]/ dist        
            vec_x = x[i_mp][0]-x[i_cell][0]
            vec_y = x[i_mp][1]-x[i_cell][1]
            vec_z = x[i_mp][2]-x[i_cell][2]       
            x_p = x[i_cell][0] + frac * vec_x
            y_p = x[i_cell][1] + frac * vec_y
            z_p = x[i_cell][2] + frac * vec_z        
            lmp.command("create_atoms 8 single " +str(x_p)+ " "+str(y_p)+ " "+str(z_p))
            n_atoms = n_atoms + 1
            dist_patch = ((x[i_cell][0]-x_p)**2.0+(x[i_cell][1]-y_p)**2.0+(x[i_cell][2]-z_p)**2.0)**0.5
            lmp.command("set atom "+str(n_atoms)+ " mol " + str(IDmp))          
            pair_cell_patch.append([ID,n_atoms,dist_patch])
            patches_IDs.append(n_atoms)
            lmp.command("create_bonds single/bond 7 " + str(ID) + " " + str(n_atoms))
        
    lmp.command("set type 8 diameter 0.1")
    lmp.command("set type 8 mass 0.000001")
    lmp.command("unfix integrator")
    lmp.command("unfix 1")
    lmp.command("unfix patches")
    lmp.command("group MPpatches delete")
    lmp.command("group langevin delete")
    lmp.command("group MPpatches type 8")
    lmp.command("group FreezeNow union FreezeNow MPpatches")
    lmp.command("group MPpatches include molecule")
    lmp.command("group langevin subtract all MPpatches zeroforce")
    lmp.command("fix integrator langevin nve")
    lmp.command("fix 1 langevin langevin ${temperature} ${temperature} ${damp} ${seed}")
    lmp.command("fix patches MPpatches rigid/small molecule langevin ${temperature} ${temperature} ${damp} ${seed}")
        
    lmp.command("run 100")
        
  #correct patches

  lmp.command("group del type 8")
  lmp.command("delete_atoms group del compress no bond yes")
  lmp.command("group del delete")    
  nlocal = lmp.extract_global("nlocal")
  x = lmp.extract_atom("x")
  ids = lmp.extract_atom("id")
  rad = lmp.extract_atom("radius")
  ids_list = []
    
  for index in range(nlocal):
    ids_list.append(ids[index])
    
  cell_patch = []
  patches_IDs = []
  distance_patchs = []
    
  n_atoms = max(ids_list)
    
  for IDmp in range(1,1201):
    i_mp = ids_list.index(IDmp)
        
    for ID in range(1201,1244):
        
      if ID in cell:
        i_cell = ids_list.index(ID) 
        Xmp = x[i_mp][0]
        X = x[i_cell][0]
        Ymp = x[i_mp][1]
        Y = x[i_cell][1]
        Zmp = x[i_mp][2]
        Z= x[i_cell][2]
        DX = (Xmp-X)**2.0
        DY = (Ymp-Y)**2.0
        DZ = (Zmp-Z)**2.0  
        dist = (DX+DY+DZ)**(0.5)
        rad_sum = rad[i_cell] + rad[i_mp]
                
        if dist <= rad_sum:
          frac = rad[i_cell]/ dist    
          vec_x = x[i_mp][0]-x[i_cell][0]
          vec_y = x[i_mp][1]-x[i_cell][1]
          vec_z = x[i_mp][2]-x[i_cell][2]      
          x_p = x[i_cell][0] + frac * vec_x
          y_p = x[i_cell][1] + frac * vec_y
          z_p = x[i_cell][2] + frac * vec_z       
          lmp.command("create_atoms 8 single " +str(x_p)+ " "+str(y_p)+ " "+str(z_p))
          n_atoms = n_atoms + 1
          dist_patch = ((x[i_cell][0]-x_p)**2.0+(x[i_cell][1]-y_p)**2.0+(x[i_cell][2]-z_p)**2.0)**0.5
          lmp.command("set atom "+str(n_atoms)+ " mol " + str(IDmp))       
          cell_patch.append(ID)
          patches_IDs.append(n_atoms)
          distance_patchs.append(dist_patch)
          lmp.command("create_bonds single/bond 7 " + str(ID) + " " + str(n_atoms))
    
  ids = lmp.extract_atom("id")
  nlocal = lmp.extract_global("nlocal")
  ids_list = []
    
  for index in range(nlocal):
    ids_list.append(ids[index])
    
  lmp.command("set type 8 diameter 0.1")
  lmp.command("set type 8 mass 0.000001")
  lmp.command("unfix integrator")
  lmp.command("unfix 1")
  lmp.command("unfix patches")
  lmp.command("group MPpatches delete")
  lmp.command("group langevin delete")
  lmp.command("group MPpatches type 8")
  lmp.command("group FreezeNow union FreezeNow MPpatches")
  lmp.command("group MPpatches include molecule")
  lmp.command("group langevin subtract all MPpatches zeroforce")
  lmp.command("fix integrator langevin nve")
  lmp.command("fix 1 langevin langevin ${temperature} ${temperature} ${damp} ${seed}")
  lmp.command("fix patches MPpatches rigid/small molecule langevin ${temperature} ${temperature} ${damp} ${seed}")
    
  lmp.command("run 100")
    
  #contraction 3

  for ID in desactive:
    if ID in layer3_ID:
      layer3_ID.remove(ID)

  distancesBegin = []
  IDBegin = []
    
  if len(layer3_ID) != 0:
    for ID in layer3_ID:
      i_c = ids_list.index(centerID)
      i_n = ids_list.index(ID)   
      centerX = x[i_c][0]
      X = x[i_n][0]
      centerY = x[i_c][1]
      Y = x[i_n][1]
      centerZ = x[i_c][2]
      Z= x[i_n][2]
      DX = (centerX-X)**2.0
      DY = (centerY-Y)**2.0
      DZ = (centerZ-Z)**2.0   
      dist = (DX+DY+DZ)**(0.5)
      distancesBegin.append(dist)     
      IDBegin.append(ID)
    
  for i in range(0,10000):
    dt = lmp.get_thermo("dt")
    dr5 = r5-r5*0.1*i*dt
    dr6 = r6-r6*0.1*i*dt
    lmp.command("bond_style harmonic")
    lmp.command("bond_coeff 1 ${bondK} "+str(dr1))
    lmp.command("bond_coeff 2 ${bondK} "+str(dr2))
    lmp.command("bond_coeff 3 ${bondK} "+str(dr3))
    lmp.command("bond_coeff 4 ${bondK} "+str(dr4))
    lmp.command("bond_coeff 5 ${bondK} "+str(dr5))
    lmp.command("bond_coeff 6 ${bondK} "+str(dr6))
    lmp.command("bond_coeff 7 ${bondmc} ${r7}") 
    lmp.command("bond_coeff 8 ${bondinit} ${r7}")
    nlocal = lmp.extract_global("nlocal")
    x = lmp.extract_atom("x")
    ids = lmp.extract_atom("id")
    ids_list = []
        
    for index in range(nlocal):
      ids_list.append(ids[index])
        
    for ID in range(1201,1244):
        
      if ID in active_layer:
            
        if ID in cell_patch:
          idx_cell = cell_patch.index(ID)
          i_cell = ids_list.index(ID)
          i_p = ids_list.index(patches_IDs[idx_cell])
          d0 = distance_patchs[idx_cell]
          cellX = x[i_cell][0]
          X = x[i_p][0]
          cellY = x[i_cell][1]
          Y = x[i_p][1]
          cellZ = x[i_cell][2]
          Z= x[i_p][2]    
          DX = (cellX-X)**2.0
          DY = (cellY-Y)**2.0
          DZ = (cellZ-Z)**2.0
          dist = (DX+DY+DZ)**(0.5)    
          diff = abs(d0-dist)
                    
          if diff >= thr_freeze:
            lmp.command("set atom " + str(ID)+" type 7")
            lmp.command("group zeroforce type 7")
            lmp.command("group langevin subtract all zeroforce MPpatches")
                        
            if ID in active_layer:
              active_layer.remove(ID)
                        
        if ID in IDBegin:
                
          if ID in active_layer:
            i_begin = IDBegin.index(ID)
            i_c = ids_list.index(centerID)
            i_n = ids_list.index(ID)       
            centerX = x[i_c][0]
            X = x[i_n][0]
            centerY = x[i_c][1]
            Y = x[i_n][1]
            centerZ = x[i_c][2]
            Z= x[i_n][2]
            DX = (centerX-X)**2.0
            DY = (centerY-Y)**2.0
            DZ = (centerZ-Z)**2.0  
            dist = (DX+DY+DZ)**(0.5)       
            thr = distancesBegin[i_begin] * 0.96
                        
            if dist < thr:
              lmp.command("set atom "+str(ID)+" type 6")
              active_layer.remove(ID)
              layer3_ID.remove(ID)
                            
              if ID not in desactive:
                desactive.append(ID)
                    
    lmp.command("group del type 6")
    lmp.command("delete_atoms group del compress no bond yes")
    lmp.command("group del delete")
        
    lmp.command("run 5")

  lmp.file("SAVElastStep.run")
    
  lmp.command("clear")

  lmp.close()
