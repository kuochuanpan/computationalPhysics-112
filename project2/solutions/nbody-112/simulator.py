import time
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from .particles import Particles
from numba import jit, njit, prange, set_num_threads

"""
The N-Body Simulator class is responsible for simulating the motion of N bodies



"""

class NBodySimulator:

    def __init__(self, particles: Particles):
        
        self.particles = particles
        self.setup() # use default settings

        return

    def setup(self, G=1,
                    rsoft=0.01,
                    method="RK4",
                    io_freq=10,
                    io_header="nbody",
                    io_screen=True,
                    visualization=False):
        """
        Customize the simulation enviroments.

        :param G: the graivtational constant
        :param rsoft: float, a soften length
        :param meothd: string, the numerical scheme
                       support "Euler", "RK2", and "RK4"

        :param io_freq: int, the frequency to outupt data.
                        io_freq <=0 for no output. 
        :param io_header: the output header
        :param io_screen: print message on screen or not.
        :param visualization: on the fly visualization or not. 
        """
        
        self.G = G
        self.rsoft = rsoft
        self.method = method.lower()
        self.io_freq = io_freq
        self.io_header = io_header
        self.io_screen = io_screen
        self.visualization = visualization


        return

    def _advance_particles(self, dt, particles):
        
        method = self.method
        if method == "euler": 
            particles = self._advance_particles_Euler(dt, particles)
        elif method == "rk2":
            particles = self._advance_particles_RK2(dt, particles)
        elif method == "rk4":
            particles = self._advance_particles_RK4(dt, particles)
        
        return particles

    def evolve(self, dt:float, tmax:float):
        """
        Start to evolve the system

        :param dt: float, the time step
        :param tmax: float, the total time to evolve
        
        """
        time = self.particles.time
        nsteps = int(np.ceil((tmax-time)/dt))
        particles = self.particles

        for n in range(nsteps):

            # make sure the last step is correct
            if (time + dt) > tmax:
                dt = tmax - time

            # updates (physics)
            particles = self._advance_particles(dt, particles)

            # check IO
            if (n % self.io_freq) == 0:
                
                # print info to screen
                if self.io_screen:
                    print("Time: ", time, "dt: ", dt)

                # check output directroy
                folder = "data_"+self.io_header
                Path(folder).mkdir(parents=True, exist_ok=True)

                # output data
                fn = self.io_header+"_"+str(n).zfill(6)+".dat"
                fn = folder+"/"+fn
                self.particles.output(fn)

                # visualization
                if self.visualization:
                    particles.draw()

            time += dt



        print("Simulation is done!")
        return

    def _calculate_acceleration(self, nparticles, masses, positions):
        """
        Calculate the acceleration of the particles
        """
        accelerations = np.zeros_like(positions)
        rsoft = self.rsoft
        G     = self.G

        # TODO
        for i in range(nparticles):
            for j in range(nparticles):
                if (j >i):
                    rij = positions[i] - positions[j]
                    r   = np.sqrt(np.sum(rij**2)+rsoft**2)
                    force = - G * masses[i,0] * masses[j,0] / r**3 * rij
                    accelerations[i] = accelerations[i] + force/masses[i,0]
                    accelerations[j] = accelerations[j] - force/masses[j,0]

        return accelerations
        
    def _advance_particles_Euler(self, dt, particles):

        nparticles = particles.nparticles
        mass = particles.masses

        pos = particles.positions
        vel = particles.velocities
        acc = self._calculate_acceleration(nparticles, mass, pos)

        pos = pos + vel*dt
        vel = vel + acc*dt
        acc = self._calculate_acceleration(nparticles, mass, pos)

        particles.set_particles(pos, vel, acc)

        return particles

    def _advance_particles_RK2(self, dt, particles):

        nparticles = particles.nparticles
        mass = particles.masses

        pos = particles.positions
        vel = particles.velocities

        t1 = time.time()
        acc = self._calculate_acceleration(nparticles, mass, pos)

        t2 = time.time()
        print("Time for acc", t2-t1)
        pos2 = pos + vel*dt
        vel2 = vel + acc*dt
        t3 = time.time()
        print("Time for updates", t3-t2)
        acc2 = self._calculate_acceleration(nparticles, mass, pos2)

        pos2 = pos2 + vel2*dt
        vel2 = vel2 + acc2*dt

        pos = 0.5*(pos + pos2)
        vel = 0.5*(vel + vel2)
        acc = self._calculate_acceleration(nparticles, mass, pos)

        particles.set_particles(pos, vel, acc)

        return particles

    def _advance_particles_RK4(self, dt, particles):
        
        #TODO








        return particles



if __name__ == "__main__":
    
    pass