import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_ivp



class SpringMassSystem:
    def __init__(self, mass, spring_constant, damping=0.0):
        self.m = mass
        self.k = spring_constant
        self.c = damping  # damping coefficient
    
    def equation_of_motion(self, t, state):
        """
        state = [position, velocity]
        Returns derivatives: [velocity, acceleration]
        """
        x, v = state
        acceleration = (-self.k * x - self.c * v) / self.m
        return [v, acceleration]
    
    def solve(self, initial_conditions, t_span, t_eval=None):
        """
        Solve the system with given initial conditions
        initial_conditions = [initial_position, initial_velocity]
        """
        solution = solve_ivp(self.equation_of_motion, t_span, 
                           initial_conditions, t_eval=t_eval)
        return solution

class AnimatedSpringMass:
    def __init__(self, mass, spring_constant, damping=0.0):
        self.system = SpringMassSystem(mass, spring_constant, damping)
        
        # Set up the figure
        self.fig, (self.ax1, self.ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Spring parameters for visualization
        self.spring_length = 2.0
        self.equilibrium_pos = 0.0
        
    def animate(self, initial_conditions, t_span, num_frames=200):
        self.t_eval = np.linspace(t_span[0], t_span[1], num_frames)
        self.solution = self.system.solve(initial_conditions, t_span, self.t_eval)
        
        # Set up animation
        self.mass_point, = self.ax1.plot([], [], 'ro', markersize=10)
        self.spring_line, = self.ax1.plot([], [], 'b-', linewidth=2)
        self.trajectory, = self.ax2.plot([], [], 'r-')
        
        # Set up axes
        self.ax1.set_xlim(-1, 1)
        self.ax1.set_ylim(-self.spring_length-1, 1)
        self.ax1.set_aspect('equal')
        self.ax1.set_title('Spring-Mass System')
        self.ax1.grid(True)
        
        self.ax2.set_xlim(t_span[0], t_span[1])
        self.ax2.set_ylim(-1.5, 1.5)
        self.ax2.set_xlabel('Time (s)')
        self.ax2.set_ylabel('Position (m)')
        self.ax2.set_title('Position vs Time')
        self.ax2.grid(True)
        
        anim = FuncAnimation(self.fig, self.update, frames=num_frames, 
                           interval=50, blit=False)
        plt.show()
        return anim
    
    def update(self, frame):
        # Current position
        x = self.solution.y[0][frame]
        
        # Update spring visualization
        spring_x = [0, 0]
        spring_y = [0, -self.spring_length - x]
        
        self.spring_line.set_data(spring_x, spring_y)
        self.mass_point.set_data([0], [-self.spring_length - x])
        
        # Update trajectory
        self.trajectory.set_data(self.t_eval[:frame+1], self.solution.y[0][:frame+1])
        
        return self.spring_line, self.mass_point, self.trajectory

# Create and run animation
animator = AnimatedSpringMass(mass=0.5, spring_constant=20.0, damping=10.0)
animation = animator.animate([1.0, 0.0], (0, 10), num_frames=200)