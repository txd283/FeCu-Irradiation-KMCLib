# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham
# This is used to convert the lattice output files to be used by the clusterCalculator.cpp
# program to calculate the data 
# . to run the code, you need to paste it at the end of the lattice
# output. A method could be 'cat convertToCluster.py >> results1.py' then 'python results1.py'
# this is a very crude way of doing this, and it loads all the data into RAM. it was quickly
# written to extract the data, if I had more time I would make it better.

#Calculates the repetitons based on the number of atoms
repetitions = (len(sites)/2)**(1/3.0)
i = 0 
# loop for the number of steps (called i).
print ("Converting the lattice trajectory to cluster.<steps>.txt files for clusterCalculator.")
for steps in range(len(steps)):
    n_atoms = len(sites)
    print("Total number of Cu  = %i"%(n_atoms*5e-4 ))
    print("Total number of Va  = 20")
    
    # open the file lattice00{}.cfg where {} is the variable which increase 0,1,2,3... by n
    steps2 = str(steps) # convert to string to be used in the filename.
    filename = 'cluster.' + steps2 + '.txt'
    file = open(filename, 'w')
    print("Writing file %s"%(filename))
    # print out the header of cfg

    j = 0
    for position in sites:
   
        s1 = position[0]/repetitions  # works
        s2 = position[1]/repetitions  # works
        s3 = position[2]/repetitions  # works
        atom_ID = types[i][j]
        
        if atom_ID == '0':
            file.write("%f %f %f\n"%(s1, s2, s3))
        elif atom_ID == '0.1':
            file.write("%f %f %f\n"%(s1, s2, s3))
        else:
            pass
        j += 1
    i += 1
    file.close()