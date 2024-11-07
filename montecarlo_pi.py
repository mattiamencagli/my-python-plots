# Try with class instead
import numpy as np
import matplotlib.pylab as plt
import matplotlib
import random 
import time
from IPython.display import display, clear_output

class MonteCarloPiEstimator:
    def __init__(self, points, dtplot):
        self.points = int(points)
        self.dtplot = int(dtplot)
        self.x = np.zeros(self.points)
        self.y = np.zeros(self.points)
        self.circle_points = 0
        self.pi_estimate = 0
        self.inside_x = []
        self.inside_y = []
        self.outside_x = []
        self.outside_y = []
        
        # Set up the figure and axis
        self.fig, self.ax = plt.subplots(figsize=(6, 6))
        self.ax.set_xlim(-1.2, 1.2)
        self.ax.set_ylim(-1.2, 1.2)
        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        self.ax.set_title(r"Estimating $\pi$ using Monte Carlo Method")
        
        # Set up the scatter plots
        self.inside_plot, = self.ax.plot([], [], 'ro', markersize=1, label='Inside Circle')
        self.outside_plot, = self.ax.plot([], [], 'bo', markersize=1, label='Outside Circle')
        
        # Draw a circle for reference
        circle = plt.Circle((0, 0), 1, color='black', fill=False, linewidth=1.5)
        self.ax.add_patch(circle)
        
        # Display the estimated pi value as text
        self.pi_text = self.ax.text(-1, 1.1, '', fontsize=12, color='purple')
        
    def run_simulation(self):
        for i in range(self.points):
            # Generate a random point
            x = random.uniform(-1, 1)
            y = random.uniform(-1, 1)
            
            # Check if the point is inside the circle
            if x**2 + y**2 <= 1:
                self.circle_points += 1
                self.inside_x.append(x)
                self.inside_y.append(y)
            else:
                self.outside_x.append(x)
                self.outside_y.append(y)
            
            if i%self.dtplot==0:
                # Update the π estimate
                pi_estimate = 4 * self.circle_points / (i + 1)
                
                # Update the scatter plot data
                self.inside_plot.set_data(self.inside_x, self.inside_y)
                self.outside_plot.set_data(self.outside_x, self.outside_y)
                
                # Update the π estimate display
                self.pi_text.set_text(f"π = {pi_estimate:.6f} - points = {i}")
                
                # Clear the previous output and display the updated plot
                # clear_output(wait=True)  # Clears the previous output to update the plot
                # display(self.fig)  # Display the updated plot
                # Explicitly draw the updated plot
                self.fig.canvas.draw()
                self.fig.canvas.flush_events() 
                plt.pause(1e-10)
            
        plt.ioff()
        plt.show()

# Example usage
estimator = MonteCarloPiEstimator(1e6, 347)
estimator.run_simulation()

