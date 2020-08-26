import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import animation
from itertools import combinations

class Particle:
    """A class representing a two-dimensional particle."""

    def __init__(self, x, y, vx, vy, radius):
        """Initialize the particle's position, velocity, and radius.

        Any key-value pairs passed in the styles dictionary will be passed
        as arguments to Matplotlib's Circle patch constructor.

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


       
     

    # For convenience, map the components of the particle's position and
    # velocity vector onto the attributes x, y, vx and vy.
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
        initialState  = "infected" if np.random.uniform(low=0.0, high=1.0) < 0.05 else "susceptible"
        if initialState == "infected":
            self.styles = {'color': 'red'}
        else:
            self.styles = {'color': 'green'}
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
        self.completeInteraction = []
        self.timeInterval = []
        self.recoveryInterval = []
        self.deathInterval = []

    def init_particles(self, n, radius):
        """Initialize the n Particles of the simulation.

        Positions and velocities are chosen randomly; radius can be a single
        value or a sequence with n values.

        """

        self.n = n
        self.particles = []
        
        while (len(self.particles) < n):
    
                # Choose x, y so that the Particle is entirely inside the
                # domain of the simulation.
                x, y = radius + (1 - 2*radius) * np.random.random(2)
                # Choose a random velocity (within some reasonable range of
                # values) for the Particle.
                vr = 0.1 * np.random.random() + 0.05
                vphi = 2*np.pi * np.random.random()
                vx, vy = vr * np.cos(vphi), vr * np.sin(vphi)
                particle = Particle(x, y, vx, vy, radius)
             
                # Check that the Particle doesn't overlap one that's already
                # been placed.
                for p2 in self.particles:
                    if p2.overlaps(particle):
                        continue
                else:
                   # print(particle.susceptibility)
                    #print(particle.covidState)
                    #print(particle.infectivity)
                    self.compartmentStats[particle.covidState] += 1
                    #print(self.compartmentStats)
                    self.particles.append(particle)

    @property
    def getCompartments(self):
        return self.compartmentStats

    def handle_collisions(self):

        """Detect and handle any collisions between the Particles.

        When two Particles collide, they do so elastically: their velocities
        change such that both energy and momentum are conserved.

        """

        def change_velocities(p1, p2):
            """
            Particles p1 and p2 have collided elastically: update their
            velocities.

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

        # We're going to need a sequence of all of the pairs of particles when
        # we are detecting collisions. combinations generates pairs of indexes
        # into the self.particles list of Particles on the fly.
        pairs = combinations(range(self.n), 2)
        for i,j in pairs:
            if self.particles[i].overlaps(self.particles[j]):

                if (self.particles[i].covidState == 'dead' or self.particles[j].covidState == 'dead'):
                    continue
                
                change_velocities(self.particles[i], self.particles[j])

                if ((self.particles[i].covidState == "infected" or self.particles[i].covidState == "asymptotic") and self.particles[j].covidState == "susceptible") :
                        temp = np.random.uniform(low=0.0, high=1.0)
                        if (temp < self.particles[i].infectivity * self.particles[j].susceptibility):
                            self.particles[j].covidState = "latent"
                            self.particles[j].timeStateStart = self.curTime
                            self.compartmentStats['latent'] += 1
                   
                            self.compartmentStats['susceptible'] -= 1
                
                if ((self.particles[j].covidState == "infected" or self.particles[j].covidState == "asymptotic") and self.particles[i].covidState == "susceptible") :
                        temp = np.random.uniform(low=0.0, high=1.0)
                        if (temp < self.particles[j].infectivity * self.particles[i].susceptibility):
                            self.particles[i].covidState = "latent"
                            self.particles[i].timeStateStart = self.curTime
                            self.compartmentStats['latent'] += 1
                            self.compartmentStats['susceptible'] -= 1

                
    @staticmethod
    def sampleRecovery(infectiontime):
        minInfectionTime = 14
        if (infectiontime < minInfectionTime):
            return 0
        else: 
            heuristic = 10.0
            return (infectiontime - minInfectionTime) / heuristic
    
    @staticmethod
    def sampleDeath(infectiontime):

        minInfectionTime = 17
        heuristic = 30.0
        if (infectiontime < 14):
            return 0
        else:
            return (infectiontime - minInfectionTime) / heuristic

    @staticmethod
    def sampleInfection(timeLatent):
        if (timeLatent < 2):
            return 0
        else:
            heuristic = 12.0
            return (timeLatent - 2.0) / heuristic
    
    @staticmethod
    def sampleAsymptotic():
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
                    p.styles = {'color': 'purple'}
                    p.covidState = "recovered"
                    p.timeStateStart = self.curTime
                
                
                elif (temp2 < self.sampleDeath(stateTimeDays) / self.ticksPerDay):
                    self.compartmentStats[p.covidState] -= 1
                    p.covidState = "dead"
                    p.styles = {'color': 'black'}
                    self.compartmentStats["dead"] += 1
                    p.timeStateStart = self.curTime
              

                
            elif (p.covidState == "latent"):
                if (temp < self.sampleInfection(stateTimeDays)/self.ticksPerDay):
                    p.timeStateStart = self.curTime
                    self.compartmentStats["latent"] -= 1
                    if (temp2 < self.sampleAsymptotic()):
                        p.covidState = "asymptomatic"
                        p.styles = {'color': 'green'}
                        self.compartmentStats["asymptomatic"] += 1
                    else:
                        p.covidState = "infected"
                        p.styles = {'color': 'red'}
                        self.compartmentStats["infected"] += 1

            if (p.covidState != "dead"):
            
                p.advance(dt)
            

            self.circles[i].center = p.r
        
        self.handle_collisions()
        return self.circles




    def init(self):
        """Initialize the Matplotlib animation."""
        self.lineDeathes.set_data([], [])
        self.lineRecovered.set_data([], [])
        self.circles = []
        for particle in self.particles:
            self.circles.append(particle.draw(self.ax))
        return self.patches


    def animate(self, i):
        """The function passed to Matplotlib's FuncAnimation routine."""
        self.curTime += 1
        self.advance_animation(0.1)
        #self.completeInteraction.append(self.compartmentStats)
        return self.circles


    def update(self, i, l1, l2):
        self.timeInterval.append(self.curTime)
        self.recoveryInterval.append(self.compartmentStats['recovered'])
        self.deathInterval.append(self.compartmentStats['dead'])

        l1.set_data(self.timeInterval, self.recoveryInterval)
        l2.set_data(self.timeInterval, self.deathInterval)
        self.curTime += 1
        self.advance_animation(1)
        return [l1,l2]

    def animateAll(self, i):
        self.timeInterval.append(self.curTime)
        self.recoveryInterval.append(self.compartmentStats['recovered'])
        self.deathInterval.append(self.compartmentStats['dead'])
        self.lineRecovered.set_data(self.timeInterval, self.recoveryInterval)
        self.lineDeathes.set_data(self.timeInterval, self.deathInterval)
        
        self.curTime += 1
        self.advance_animation(0.1)

        return self.patches

    

    def do_animation(self, save=False):
        """Set up and carry out the animation of the molecular dynamics.

        To save the animation as a MP4 movie, set save=True.
        """

        
        fig =  plt.figure()
        self.ax = fig.add_subplot(2,2,1)
        for s in ['top','bottom','left','right']:
            self.ax.spines[s].set_linewidth(2)
        self.ax.set_aspect('equal', 'box')
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)
        self.ax.xaxis.set_ticks([])
        self.ax.yaxis.set_ticks([])
        self.ax_line = fig.add_subplot(1,2,2)
        self.ax_line.set_ylim([0,100])
        self.ax_line.set_xlim([0,500])
     
        
        self.lineRecovered = self.ax_line.plot([],[], color ="r")[0]
        self.lineDeathes = self.ax_line.plot([],[], color="g")[0]

        self.patches = [self.ax, self.lineRecovered, self.lineDeathes]

        anim = animation.FuncAnimation(fig, self.animateAll, init_func=self.init,
                               frames=100, interval=20, blit=True)





        """  self.ax = fig.add_subplot(1,2,2)
        self.ax_line = fig.add_subplot(2,2,1)
        l1, = self.ax_line.plot([], [], color = "r")
        l2, = self.ax_line.plot([], [], color = "g")


        self.ax_line.set_ylim([0,100])
        self.ax_line.set_xlim([0,500])
        for s in ['top','bottom','left','right']:
            self.ax.spines[s].set_linewidth(2)
        self.ax.set_aspect('equal', 'box')
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)
        self.ax.xaxis.set_ticks([])
        self.ax.yaxis.set_ticks([])
        anim = animation.FuncAnimation(fig, self.animate, init_func=self.init,
                               frames=100, interval=2, blit=True)
        ani = animation.FuncAnimation(fig, self.update, 100, fargs=[l1,l2],
                  interval=2, blit=True) """
        
        if save:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=100, bitrate=1800)
            anim.save('collision.mp4', writer=writer)
        else:
            plt.show()
                    


if __name__ == '__main__':
    nparticles = 100
    #radii = np.random.random(nparticles)*0.03+0.02
    sim = World(nparticles, 0.01, 10)
    sim.do_animation(save=False)