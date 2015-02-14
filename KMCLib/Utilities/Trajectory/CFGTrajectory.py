""" Module for the XYZTrajectory object """


# Copyright (c)  2013  Mikael Leetmaa
#
# This file is part of the KMCLib project distributed under the terms of the
# GNU General Public License version 3, see <http://www.gnu.org/licenses/>.
#


import sys
import time


from KMCLib.Backend.Backend import MPICommons
from KMCLib.Utilities.Trajectory.Trajectory import Trajectory


class CFGTrajectory(Trajectory):
    
    #Class for handling cfg IO to a trajectory file.
    
    def __init__(self,
                 trajectory_filename,
                 configuration,
                 max_buffer_size=None,
                 max_buffer_time=None):
        """
        Constructor for the CFGTrajectory.

        :param trajectory_filename: The file name to write trajectory information to.
        :type trajectory_filename: str

        :param sites: The lattice sites of the configuration as an Nx3 list.

        :param max_buffer_size: The max size of the the buffer in memory
                                before writing to file.
        :type max_buffer_size: int

        :param max_buffer_time: The max time limit between dumps to file.
        :type max_buffer_time: float
        """
        # Call the base class constructor.
        Trajectory.__init__(self,
                            trajectory_filename,
                            max_buffer_size,
                            max_buffer_time)

        # Init the member data.
        self.__atom_id_types = []
        self.__atom_id_coordinates = []
        self.__time = []
        self.__step = []

    def _storeData(self, simulation_time, step, configuration):
        """
        Append the coordinate, types and time information to the
        internal buffers.

        :param simulation_time: The current time of the simulation.
        :type simulation_time: float

        :param step: The step number in the simulation.
        :type step: int

        :param configuration: The configuration of the simulation.
        """
        # Store to the local buffers.
        self.__time.append(simulation_time)
        self.__step.append(step)
        self.__atom_id_coordinates.append(configuration.atomIDCoordinates())
        self.__atom_id_types.append(configuration.atomIDTypes())

    def _bufferSize(self):
        """
        Calculate and return the buffer size.
        """
        size =  sys.getsizeof(self.__atom_id_coordinates)
        size +=  sys.getsizeof(self.__atom_id_coordinates[0])*len(self.__atom_id_coordinates)
        size += sys.getsizeof(self.__atom_id_types)
        size += sys.getsizeof(self.__atom_id_types[0])*len(self.__atom_id_types)
        size += sys.getsizeof(self.__time)
        size += sys.getsizeof(self.__time[0])*len(self.__time)
        size += sys.getsizeof(self.__step)
        size += sys.getsizeof(self.__step[0])*len(self.__step)

        return size

    def flush(self):
        """ Write all buffers to file. """
        if not len(self.__step) < 1:

            # Make sure only master writes.
            if MPICommons.isMaster():

                # Write data to file.
                with open(self._trajectory_filename, 'a') as trajectory:
                    for i in range(len(self.__step)):

                        step = self.__step[i]
                        time = self.__time[i]
                        n_atoms = len(self.__atom_id_types[i])
                        
                        # calculate the super-cell length
                        repetitions = 4.0
                        H = 2.87 * repetitions
                        
                        #supercell headers for CFG file for atomeye
                        trajectory.write("Number of particles = 128\n\n")
                        trajectory.write("A = 10.0 Angstrom\n")
                        trajectory.write("H0(1,1) = %f A\n"%(H))
                        trajectory.write("H0(1,2) = 0 . 0 A\n")  
                        trajectory.write("H0(1,3) = 0 . 0 A\n") 
                        trajectory.write("H0(2,1) = 0 . 0 A\n") 
                        trajectory.write("H0(2,2) = %f A\n"%(H))
                        trajectory.write("H0(2,3) = 0 . 0 A\n")  
                        trajectory.write("H0(3,1) = 0 . 0 A\n") 
                        trajectory.write("H0(3,2) = 0 . 0 A\n")
                        trajectory.write("H0(3,3) = %f A\n"%(H))
                        trajectory.write("#\n")
                                
                        for j in range(n_atoms):
                            t = self.__atom_id_types[i][j]
                            c = self.__atom_id_coordinates[i][j]
                            
                            #reduce the coordinates
                            d = c[0]/H
                            e = c[1]/H
                            g = c[2]/H
                            
                            #change the mass atomic number depending on type of atoms there
                            if t == 'Fe':
                                mass = 55.845
                            elif t == 'V':
                                mass = 1.0
                            
                            # write the cfg file with reduced coordinates
                            trajectory.write("%f %s %f %f %f 0.0 0.0 0.0\n"%(mass, t, d, e, g))
                        
                        
            # While the other processes wait.
            MPICommons.barrier()

            # Reset the buffers.
            self.__atom_id_types = []
            self.__atom_id_coordinates = []
            self.__time = []
            self.__step = []
