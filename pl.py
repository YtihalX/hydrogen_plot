import matplotlib.pyplot as plt
import numpy as np

# Load the data from 'wavefunction' file
data = np.loadtxt('wavefunction')
x, y, z, intensity = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

# Create the plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Scatter plot without edges
sc = ax.scatter(x, y, z, c=intensity, cmap='cividis', edgecolors='none', s = 0.2)

# Add colorbar
cbar = plt.colorbar(sc)

# Labels and titles
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
plt.title('hydrogen')
ax.view_init(elev=30)
# Save the figure
plt.savefig("output_plot_python.png", dpi = 330)

# Show the plot
# plt.show()

