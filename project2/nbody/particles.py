import numpy as np
import matplotlib.pyplot as plt


class Particles:
    """
    Particle class to store particle properties
    """
    def __init__(self, N:int = 100):
        """
        :param N: number of particles.
        """
        self.nparticles     = N
        self._masses        = np.ones((N,1))
        self._positions     = np.zeros((N,3))
        self._velocities    = np.zeros((N,3))
        self._accelerations = np.zeros((N,3))
        self._tags          = np.arange(N)
        self.time = 0
    pass

    @property
    def tags(self):
        return self._tags
    
    @tags.setter
    def tags(self,some_name):
        if len(some_name) != self.nparticles:
            print("Name is too long or too short!!")
            raise ValueError
        self._tags = some_name
    
    @property
    def masses(self):
        return self._masses
    
    @masses.setter
    def masses(self,some_masses):
        size = some_masses.shape
        if size[1] != 1:
            print("Mass should be a N*1 array.")
            raise ValueError
        self._tags = some_masses
        return
    
    @property
    def positions(self):
        return self._positions
    
    @positions.setter
    def positions(self,some_x):
        if len(some_x) != self.nparticles or some_x.shape[1] != 3:
            print("positions is too long or too short!!")
            raise ValueError
        self._positions = some_x

    @property
    def velocities(self):
        return self._velocities
    
    @velocities.setter
    def velocities(self,some_v):
        if len(some_v) != self.nparticles or some_v.shape[1] != 3:
            print("velocities is too long or too short!!")
            raise ValueError
        self._velocities = some_v
        return

    @property
    def accelerations(self):
        return self._accelerations
    
    @accelerations.setter
    def accelerations(self,some_a):
        if len(some_a) != self.nparticles or some_a.shape[1] != 3:
            print("accelerations is too long or too short!!")
            raise ValueError
        self._accelerations = some_a

    def add_particles(self,masses, positions, velocities, accelerations):
        M      = self._masses
        X      = self._positions
        V      = self._velocities
        A      = self._accelerations

        number = len(masses)
        self.nparticles += number
        self._masses = np.concatenate((M,masses),axis=0)
        self._positions = np.concatenate((X , positions),axis=0)
        self._velocities = np.concatenate((V, velocities),axis=0)
        self._accelerations = np.concatenate((A, accelerations),axis=0)
        return
    
    def output(self,filename):
        f = open(filename,'w')
        f.write("#  {}  {}  {}  {}  {} \n".format("tags","masses","positions","velocities","acceleration"))
        f.write("#  {}  {}  {}  {}  {} \n".format(self.tags,self.masses,self.positions,self.velocities,self.accelerations))


if __name__ =='__main__':
    num_particles = 100
    masses        = np.ones((num_particles,1))
    print(len(masses))