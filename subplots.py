import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import animation
from itertools import combinations
import random

class Particle:
    """A generic person with attributes in the world simulation"""

    def __init__(self, x, y, vx, vy, radius):
        """Initialize the particle's pos, vel, and susceptiability and state factors

        """

        self.r = np.array((x, y))
        self.v = np.array((vx, vy))
        self.radius = radius
        self.styles = {}
        self.covidState = self.sampleCovidState()
        self.susceptibility = self.sampleSusceptibility()
        self.infectivity = self.sampleInfectivity()
        self.age = 20
        self.timeStateStart = 0



    @property
    def x(self):
        return self.r[0]
    @x.setter
    def x(self, value):
        self.r[0] = value
    @property
    def y(self):
        return self.r[1]
    @y.setter
    def y(self, value):
        self.r[1] = value
    @property
    def vx(self):
        return self.v[0]
    @vx.setter
    def vx(self, value):
        self.v[0] = value
    @property
    def vy(self):
        return self.v[1]
    @vy.setter
    def vy(self, value):
        self.v[1] = value


    

    @staticmethod
    def randn_bm(mean, std):
        u = 0
        v = 0
        while(u == 0): 
            u = np.random.uniform(low=0.0, high=1.0)
        while(v == 0): 
            v = np.random.uniform(low=0.0, high=1.0)
        pre_transform = np.sqrt( -2.0 * np.log( u ) ) * np.cos( 2.0 * np.pi * v )
        return pre_transform * std + mean
    
    
    def sampleCovidState(self):   
        initialState  = "infected" if random.random() < 0.05 else "susceptible"
        if initialState == "infected":
            self.styles = {'color': '#F07C4E'}
        else:
            self.styles = {'color': '#B0C4DE'}
        return initialState
    
    def sampleSusceptibility(self):
        return self.randn_bm(1, .1)
    
    def sampleInfectivity(self):
        return self.randn_bm(1, .1)

    def overlaps(self, other):
        """Does the circle of this Particle overlap that of other?"""

        return np.hypot(*(self.r - other.r)) < self.radius + other.radius


    def draw(self, ax):
        """Add this Particle's Circle patch to the Matplotlib Axes ax."""

        circle = Circle(xy=self.r, radius=self.radius, **self.styles)
        ax.add_patch(circle)
        return circle

    def advance(self, dt):

        """Advance the Particle's position forward in time by dt."""

        self.r += self.v * dt

        # Make the Particles bounce off the walls
        if self.x - self.radius < 0:
            self.x = self.radius
            self.vx = -self.vx
        if self.x + self.radius > 1:
            self.x = 1-self.radius
            self.vx = -self.vx
        if self.y - self.radius < 0:
            self.y = self.radius
            self.vy = -self.vy
        if self.y + self.radius > 1:
            self.y = 1-self.radius
            self.vy = -self.vy


class World:

    def __init__(self, n, radius, ticksPerDay):     
        self.compartmentStats = {"latent":0, "asymptomatic":0, "infected":0, "susceptible": 0, "recovered":0, "dead":0}
        self.init_particles(n, radius)
        
        self.curTime = 0
        self.ticksPerDay = ticksPerDay
        self.timeInterval = []
        self.recoveryInterval = []
        self.deathInterval = []
        self.asymptomaticInterval = []
        self.susceptibleInterval = []
        self.latentInterval = []
        self.infectedInterval = []

    def init_particles(self, n, radius):
        """Intialize particles from particle class with random pos, vel, 
        and a state

        """

        self.n = n
        self.particles = []
        
        while (len(self.particles) < n):
    
        
                x, y = radius + (1 - 2*radius) * np.random.random(2)
               
        
                vr = 0.1 * np.random.random() + 0.05
                vphi = 2*np.pi * np.random.random()
                vx, vy = vr * np.cos(vphi), vr * np.sin(vphi)
              
                particle = Particle(x, y, vx, vy, radius)
             
                """Removes overlapping particles """
                for p2 in self.particles:
                    if p2.overlaps(particle):
                        continue
                else:
                 
                    self.compartmentStats[particle.covidState] += 1
                    
                    self.particles.append(particle)




    def handle_collisions(self):

        """Function to handle collisions and allocate a state based on covidState of interacting
        particles

        """

        def change_velocities(p1, p2):
            """
            Use momentum and velocity to pick decide velocity of two particles after ellastic collision

            """

            m1, m2 = p1.radius**2, p2.radius**2
            M = m1 + m2
            r1, r2 = p1.r, p2.r
            d = np.linalg.norm(r1 - r2)**2
            v1, v2 = p1.v, p2.v
            u1 = v1 - 2*m2 / M * np.dot(v1-v2, r1-r2) / d * (r1 - r2)
            u2 = v2 - 2*m1 / M * np.dot(v2-v1, r2-r1) / d * (r2 - r1)
            p1.v = u1
            p2.v = u2

        """ Create pairs to evaluate whether any two overlapped """ 
        pairs = combinations(range(self.n), 2)
        for i,j in pairs:
            if self.particles[i].overlaps(self.particles[j]):

                if (self.particles[i].covidState == 'dead' or self.particles[j].covidState == 'dead'):
                    continue
                
                change_velocities(self.particles[i], self.particles[j])

                if ((self.particles[i].covidState == "infected" or self.particles[i].covidState == "asymptomatic") and self.particles[j].covidState == "susceptible") :
                        temp = random.random()
                        if (temp < self.particles[i].infectivity * self.particles[j].susceptibility):
                            self.particles[j].covidState = "latent"
                            self.particles[j].styles = {'color': '#FAD9B9'}
                            self.particles[j].timeStateStart = self.curTime
                            self.compartmentStats['latent'] += 1  
                            self.compartmentStats['susceptible'] -= 1
                
                if ((self.particles[j].covidState == "infected" or self.particles[j].covidState == "asymptomatic") and self.particles[i].covidState == "susceptible") :
                        temp = random.random()
                        if (temp < self.particles[j].infectivity * self.particles[i].susceptibility):

                            self.particles[i].covidState = "latent"
                            self.particles[i].styles = {'color': '#FAD9B9'}
                            self.particles[i].timeStateStart = self.curTime
                            self.compartmentStats['latent'] += 1
                            self.compartmentStats['susceptible'] -= 1

                
    @staticmethod
    def sampleRecovery(infectiontime):
        """ Returns probability of recovery at certain point of time in infection state """
        minInfectionTime = 14
        if (infectiontime < minInfectionTime):
            return 0
        else: 
            heuristic = 10.0
            return (infectiontime - minInfectionTime) / heuristic
    
    @staticmethod
    def sampleDeath(infectiontime):
        """ Returns porbability of death with input as time infected """
        minInfectionTime = 17
        heuristic = 30.0
        if (infectiontime < 14):
            return 0
        else:
            return (infectiontime - minInfectionTime) / heuristic

    @staticmethod
    def sampleInfection(timeLatent):
        """ Returns infection probability and requires time patient was latent as input
        """ 
        if (timeLatent < 2):
            return 0
        else:
            heuristic = 12.0
            return (timeLatent - 2.0) / heuristic
    
    @staticmethod

    def sampleAsymptotic():
        """ Returns probability of being asymptomatic """
        return 0.1                    
        
    def advance_animation(self, dt):
        for i, p in enumerate(self.particles):
            temp = np.random.uniform(low=0.0, high=1.0)
            temp2 = np.random.uniform(low=0.0, high=1.0)

            stateTimeDays = (self.curTime - p.timeStateStart) / self.ticksPerDay
       
            if p.covidState == "infected" or p.covidState == "asymptomatic":
                if (temp < self.sampleRecovery(stateTimeDays)/self.ticksPerDay):
                    self.compartmentStats[p.covidState] -= 1
                    self.compartmentStats["recovered"] += 1
                    p.styles = {'color': '#9370DB'}
                    p.covidState = "recovered"
                    p.timeStateStart = self.curTime
                
                
                elif (temp2 < self.sampleDeath(stateTimeDays) / self.ticksPerDay):
                    self.compartmentStats[p.covidState] -= 1
                    p.covidState = "dead"
                    p.styles = {'color': '#000000'}
                    self.compartmentStats["dead"] += 1
                    p.timeStateStart = self.curTime
              

                
            elif (p.covidState == "latent"):
                if (temp < self.sampleInfection(stateTimeDays)/self.ticksPerDay):
                    p.timeStateStart = self.curTime
                    self.compartmentStats["latent"] -= 1
                    if (temp2 < self.sampleAsymptotic()):
                        p.covidState = "asymptomatic"
                        p.styles = {'color': '#F39F79'}
                        self.compartmentStats["asymptomatic"] += 1
                    else:
                        p.covidState = "infected"
                        p.styles = {'color': '#F07C4E'}
                        self.compartmentStats["infected"] += 1

            if (p.covidState != "dead"):
            
                p.advance(dt)
            

            self.circles[i].center = p.r
        
        self.handle_collisions()
        return self.circles




    def init(self):
        """ Initialize the lines and circles with empty arrays and initial positions 
        respectively
        """
        self.lineDeathes.set_data([], [])
        self.lineRecovered.set_data([], [])
        self.lineAsymptotic.set_data([], [])
        self.lineSusceptible.set_data([],[])
        self.lineLatent.set_data([],[])
        self.lineInfected.set_data([],[])
        self.circles = []
        for particle in self.particles:
            self.circles.append(particle.draw(self.ax))
        return self.patches


    def animate(self, i):
        """Animation for just moving the particles forward"""
        self.curTime += 1
        self.advance_animation(0.1)
     
        return self.circles

    def animateAll(self, i):
        
        self.timeInterval.append(self.curTime)
        self.recoveryInterval.append(self.compartmentStats['recovered'])
        self.deathInterval.append(self.compartmentStats['dead'])
        self.asymptomaticInterval.append(self.compartmentStats['asymptomatic'])
        self.susceptibleInterval.append(self.compartmentStats['susceptible'])
        self.latentInterval.append(self.compartmentStats['latent'])
        self.infectedInterval.append(self.compartmentStats['infected'])
        
        self.lineRecovered.set_data(self.timeInterval, self.recoveryInterval)
        self.lineDeathes.set_data(self.timeInterval, self.deathInterval)
        self.lineSusceptible.set_data(self.timeInterval, self.susceptibleInterval)
        self.lineLatent.set_data(self.timeInterval, self.latentInterval)
        self.lineAsymptotic.set_data(self.timeInterval, self.asymptomaticInterval)
        self.lineInfected.set_data(self.timeInterval, self.infectedInterval)
        
        self.ax_line.set_xlim(0, max(self.timeInterval))
        
        
        self.curTime += 1
        self.advance_animation(0.1)

        return self.patches

    

    def do_animation(self, save=False):
        """Set up lines subplot and animation for particles side-by-side and give
        saving options
        """

        
        fig =  plt.figure()
        plt.tight_layout()
        self.ax = fig.add_subplot(221)
        for s in ['top','bottom','left','right']:
            self.ax.spines[s].set_linewidth(2)
        self.ax.set_aspect('equal', 'box')
       
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)
        self.ax.xaxis.set_ticks([])
        self.ax.yaxis.set_ticks([])
        self.ax_line = fig.add_subplot(222)
        self.ax_line.set_ylim(0, self.n)

        self.lineSusceptible = self.ax_line.plot([],[], color="#B0C4DE", label="susceptible")[0]
        self.lineLatent = self.ax_line.plot([],[], color="#FAD9B9", label="latent")[0]
        self.lineAsymptotic = self.ax_line.plot([],[], color="#F39F79", label="asympotmatic")[0]
        self.lineInfected = self.ax_line.plot([],[], color="#F07C4E" ,label="infected")[0]
        self.lineRecovered = self.ax_line.plot([],[], color ="#9370DB", label="recovered")[0]
        self.lineDeathes = self.ax_line.plot([],[], color="#000000", label="dead")[0]

        self.patches = [self.ax, self.lineSusceptible, self.lineLatent, self.lineAsymptotic, self.lineRecovered, self.lineDeathes]
        
        anim = animation.FuncAnimation(fig, self.animateAll, init_func=self.init,
                               frames=100, interval=20, blit=False)
        
        if save:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=100, bitrate=1800)
            anim.save('CovidModel.mp4', writer=writer)
        else:
            plt.legend()
            plt.show()
                    


if __name__ == '__main__':
    nparticles = 100
    sim = World(nparticles, 0.01, 10)
    sim.do_animation(save=False)