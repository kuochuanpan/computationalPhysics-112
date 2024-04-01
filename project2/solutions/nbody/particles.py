import numpy as np
import matplotlib.pyplot as plt

"""

The Particles class handle all particle properties

    for the N-body simulation. 

    
Author: Kuo-Chuan Pan, NTHU 2022.10.30
                            2024.03.31
For the course, computational physics, NTHU
"""

class Particles:
    """
    Particle class to store particle properties
    """
    def __init__(self, N:int=0):
        """
        Allocate memories for particle properties

        :param N: number of particles
        """
        self.nparticles = N 
        self._masses = np.ones((N,1))
        self._positions = np.zeros((N,3))
        self._velocities = np.zeros((N,3))
        self._accelerations = np.zeros((N,3))
        self._tags = np.arange(N)
        self._time = 0.0
        return
    
    @property
    def masses(self):
        return self._masses
    
    @property
    def positions(self):
        return self._positions
    
    @property
    def velocities(self):
        return self._velocities
    
    @property
    def accelerations(self):
        return self._accelerations
    
    @property
    def tags(self):
        return self._tags
    
    @property
    def time(self):
        return self._time
    
    @masses.setter
    def masses(self, m: np.ndarray):
        if m.shape != np.ones((self.nparticles,1)).shape:
            print("Number of particles does not match!")
            raise ValueError
        
        self._masses = m
        return
    
    @positions.setter
    def positions(self, pos: np.ndarray):
        if pos.shape != np.zeros((self.nparticles,3)).shape:
            print("Number of particles does not match!")
            raise ValueError
        
        self._positions = pos
        return
    
    @velocities.setter
    def velocities(self, vel: np.ndarray):
        if vel.shape != np.zeros((self.nparticles,3)).shape:
            print("Number of particles does not match!")
            raise ValueError
        
        self._velocities = vel
        return
    
    @accelerations.setter
    def accelerations(self, acc: np.ndarray):
        if acc.shape != np.zeros((self.nparticles,3)).shape:
            print("Number of particles does not match!")
            raise ValueError
        
        self._accelerations = acc
        return
    
    @tags.setter
    def tags(self, tag: np.ndarray):
        if tag.shape != np.arange(self.nparticles).shape:
            print("Number of particles does not match!")
            raise ValueError
        
        self._tags = tag
        return

    def set_particles(self, pos, vel, acc):
        """
        Set particle properties for the N-body simulation

        :param pos: positions of particles
        :param vel: velocities of particles
        :param acc: accelerations of particles
        """
        self.positions = pos
        self.velocities = vel
        self.accelerations = acc
        return
    
    def add_particles(self, mass, pos, vel, acc):
        """
        Add N particles to the N-body simulation at once

        :param pos: positions of particles
        :param vel: velocities of particles
        :param acc: accelerations of particles
        """
        self.nparticles += mass.shape[0]
        self.masses = np.vstack((self.masses, mass))
        self.positions = np.vstack((self.positions, pos))
        self.velocities = np.vstack((self.velocities, vel))
        self.accelerations = np.vstack((self.accelerations, acc))
        self.tags = np.arange(self.nparticles)
        return

    def output(self, filename):
        """
        Output particle properties to a file

        :param filename: output file name
        """
        masses = self.masses
        pos = self.positions
        vel = self.velocities
        acc = self.accelerations
        tags = self.tags
        time = self.time

        header = """
                ----------------------------------------------------
                Data from a 3D direct N-body simulation. 

                rows are i-particle; 
                coumns are :mass, tag, x ,y, z, vx, vy, vz, ax, ay, az

                NTHU, Computational Physics 

                ----------------------------------------------------
                """
        header += "Time = {}".format(time)
        np.savetxt(filename,(tags[:],masses[:,0],pos[:,0],pos[:,1],pos[:,2],
                            vel[:,0],vel[:,1],vel[:,2],
                            acc[:,0],acc[:,1],acc[:,2]),header=header)
        return
    
    def draw(self, dim=2):
        """
        Draw particles in 3D space
        """
        fig = plt.figure()

        if dim == 2:
            ax = fig.add_subplot(111)
            ax.scatter(self.positions[:,0], self.positions[:,1])
            ax.set_xlabel('X [code unit]')
            ax.set_ylabel('Y [code unit]')
            
        elif dim == 3:
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(self.positions[:,0], self.positions[:,1], self.positions[:,2])
            ax.set_xlabel('X [code unit]')
            ax.set_ylabel('Y [code unit]')
            ax.set_zlabel('Z [code unit]')
        else:
            print("Invalid dimension!")
            return

        ax.set_aspect('equal')
        plt.tight_layout()
        plt.show()
        return fig, ax
    
    


