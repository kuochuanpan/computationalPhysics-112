import numpy as np
import matplotlib.pyplot as plt


class Particles:
    """
    Particle class to store particle properties
    """
    def __init__(self, N:100):
        self.nparticles = N
        self.time = 0
        self._masses = np.ones((N,1))
        self._positions = np.zeros((N,3))
        self._velocities = np.zeros((N,3))
        self._accelerations = np.zeros((N,3))
        self._tags = np.linspace(1,N,N)
        return
    
    @property
    def masses(self):
        return self._masses
    
    @masses.setter
    def masses(self, masses):
        if masses.shape != (self.nparticles,1):
            raise ValueError("Masses must be a vector of size N")
        self._masses = masses
        return
    
    @property
    def positions(self):
        return self._positions
    
    @positions.setter
    def positions(self, positions):
        if positions.shape != (self.nparticles,3):
            raise ValueError("Positions must be a Nx3 matrix")
        self._positions = positions
        return
    
    @property
    def velocities(self):
        return self._velocities
    
    @velocities.setter
    def velocities(self, velocities):
        if velocities.shape != (self.nparticles,3):
            raise ValueError("Velocities must be a Nx3 matrix")
        self._velocities = velocities
        return
    
    @property
    def accelerations(self):
        return self._accelerations
    
    @accelerations.setter
    def accelerations(self, accelerations):
        if accelerations.shape != (self.nparticles,3):
            raise ValueError("Accelerations must be a Nx3 matrix")
        self._accelerations = accelerations
        return

    @property
    def tags(self):
        return self._tags

    @tags.setter
    def tags(self, tags):
        if len(tags) != self.nparticles:
            raise ValueError("Number of tags must be equal to number of particles")
        self._tags = tags
        return
    
    def output(self, filename):
        """
        Output particle properties to a text file
        """
        header = "# time,tag,mass,x,y,z,vx,vy,vz,ax,ay,az\n"
        header+= "# s,,kg,m,m,m,m/s,m/s,m/s,m/s^2,m/s^2,m/s^2\n"

        np.savetxt(filename, 
                   np.hstack((np.ones((self.nparticles,1))*self.time, 
                                       self._tags.reshape(-1,1), 
                                       self._masses, 
                                       self._positions, 
                                       self._velocities, 
                                       self._accelerations)), 
                                       delimiter=",", header=header, comments="")

        return
    
    def draw(self, dim=2, save=False, filename="particles.png"):
        """
        Draw particle positions
        """
        if dim == 2:
            fig = plt.figure()
            plt.scatter(self.positions[:,0], self.positions[:,1])
            plt.xlabel("x")
            plt.ylabel("y")
            # set aspect ratio to be equal
            plt.gca().set_aspect('equal', adjustable='box')
            plt.tight_layout()
            plt.show()
            return fig, ax
        elif dim == 3:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(self.positions[:,0], self.positions[:,1], self.positions[:,2])
            # set aspect ratio to be equal
            ax.set_aspect('equal')
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("z")
            plt.tight_layout()
            plt.show()
            return fig, ax
        else:
            raise ValueError("Invalid dimension")
        
    
    
if __name__=='__main__':

    pts = Particles(100)
    tags = pts.tags
    pts.tags = [-1]
    print(tags)


