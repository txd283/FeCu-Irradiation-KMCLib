
# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham
# This is used to convert the lattice output files to CFG files for atomeye
# visualistion. to run the code, you need to paste it at the end of the lattice
# output. A method could be 'cat convertToCFG.py >> results1.py' then 'python results1.py'
# this is a very crude way of doing this, and it loads all the data into RAM. it was quickly
# written to extract the data, if I had more time I would make it better.

# Calculates the repetitons based on the number of atoms
repetitions = (len(sites)/2)**(1/3.0)
mass = 0.0
i = 0
 
# loop for the number of steps (called i).
print ("Converting the lattice trajectory to .cfg files for atomeye.")

for steps in range(len(steps)):
    n_atoms = len(sites)
    print("Total number of atoms = %i"%(n_atoms))
    print ("Step = %f"%(steps))
    
    # open the file lattice00{}.cfg where {} is the variable which increase 0,1,2,3... by n
    steps2 = str(steps) # convert to string to be used in the filename.
    filename = 'cfg0' + steps2 + '.cfg'
    file = open(filename, 'w')
    print("Writing file %s"%(filename))
    
    # print out the header of cfg
    file.write("Number of particles = %i\n\n"%(n_atoms))
    file.write("A = 4 Angstrom\n")
    file.write("H0(1,1) = %.1f A\n"%(repetitions))
    file.write("H0(1,2) = 0 . 0 A\n")
    file.write("H0(1,3) = 0 . 0 A\n")
    file.write("H0(2,1) = 0 . 0 A\n")
    file.write("H0(2,2) = %.1f A\n"%(repetitions))
    file.write("H0(2,3) = 0 . 0 A\n")
    file.write("H0(3,1) = 0 . 0 A\n")
    file.write("H0(3,2) = 0 . 0 A\n")
    file.write("H0(3,3) = %.1f A\n"%(repetitions))
    file.write("#\n")
    file.write("# step = %.0f\n"%(steps))
    
    # loop the types of atoms array
    j = 0
    for position in sites:
    
        s1 = position[0]/repetitions  # works
        s2 = position[1]/repetitions  # works
        s3 = position[2]/repetitions  # works
        atom_ID = types[i][j]
        if atom_ID == '1':
            mass = 55.845
            type_atom = 'Fe'
        elif atom_ID == '0':
            mass = 1.0
            type_atom = 'H'
        elif atom_ID == '0.1':
            mass = 63.543
            type_atom = 'Cu'
        file.write("%.2f %s %.2f %.2f %.2f 0.0 0.0 0.0\n"%(mass, type_atom, s1, s2, s3))
        j += 1
    i += 1
    file.close()