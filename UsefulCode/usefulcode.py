# Author = Thomas Davis, email = txd283@bham.ac.uk / University of Birmingham
# Used to calculate the nearest neighbours in a supercell of 2.0 cutoff in CustomRate
"""
  NN1 = numpy.zeros(shape=(15, 8))
  NN2 = numpy.zeros(shape=(15, 6))
  for j in range(0,15):
      nn1count = 0
      nn2count = 0
      for k in range(0,65):
          #if k == j:
          #    pass
          distance_sq = (geometry[j][0] - geometry[k][0])**2 + (geometry[j][1] - geometry[k][1])**2 + (geometry[j][2] - geometry[k][2])**2
          if distance_sq == 0.75:
              NN1[j][nn1count] = k
              nn1count += 1
          elif distance_sq == 1.0:
              NN2[j][nn2count] = k
              nn2count += 1
          #else:
          #    pass
          
  print(NN1)
  print(NN2)	   
"""
	   
	   