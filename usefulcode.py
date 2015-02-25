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
	   
	   
	   
       b_first_neighbours_Fe = int(types_before[1]) + int(types_before[2]) + int(types_before[3]) + int(types_before[4]) + int(types_before[5]) + int(types_before[6]) + int(types_before[7]) + int(types_before[8])
       
       b_first_neighbours_Va = 8 - b_first_neighbours_Fe
       
       a_first_neighbours_Fe = int(types_after[1]) + int(types_after[2]) + int(types_after[3]) + int(types_after[4]) + int(types_after[5]) + int(types_after[6]) + int(types_after[7]) + int(types_after[8])
       a_first_neighbours_Va = 8 - a_first_neighbours_Fe
       
       b_second_neighbours_Fe = int(types_before[9]) + int(types_before[10]) + int(types_before[11]) + int(types_before[12]) + int(types_before[13]) + int(types_before[14])
       b_second_neighbours_Va = 6 - b_second_neighbours_Fe
       
       a_second_neighbours_Fe = int(types_after[9]) + int(types_after[10]) + int(types_after[11]) + int(types_after[12]) + int(types_after[13]) + int(types_after[14])
       a_second_neighbours_Va = 6 - a_second_neighbours_Fe
      
	   
	   
	   