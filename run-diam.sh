#!/bin/bash

python3 cell_disk_data.py

#diameter of the microparticle
size=1.00

#for each bed file
for bed in {1..16}

do 

  #for cell file $1 to $2
  for ((k=$1; k<=$2; k++))

  do

    ########## MD part ###############

    #removes the old 0.model.py
    rm -f 0.model.py
		
    #replaces in model.py ID_FILE SIZE and BEDID and saves in a new file 0.model.py
    cat model.py | sed 's/ID_FILE/'${k}'/;s/SIZE/'${size}'/;s/BEDID/'${bed}'/'>> 0.model.py

    ###### submition part ############
		
    #run 0.model.py in python 2.7
    python 0.model.py	

  done 

done
