using Plots
using DelimitedFiles

plotlyjs()

# Load data from a file, assuming it's space-delimited
data = readdlm("wavefunction")
x, y, z, prob = data[:,1], data[:,2], data[:,3], data[:,4]

# Create the plot
p = scatter3d(x, y, z, colorbar=true, marker_z=prob,
              xlabel="X-axis", ylabel="Y-axis", zlabel="Z-axis", title="4D Plot", dpi = 300, markersize = 1, markerstrokestyle = :dot, clims = (-2, 10))

# Save the plot
savefig(p, "output_plot.png")

