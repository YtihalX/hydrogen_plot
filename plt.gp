# Set the terminal
set style fill transparent solid 0.5

set terminal pngcairo size 800,800

# Set the output file
set output 'output1.png'
set xlabel 'x'
set ylabel 'y'
set zlabel 'z'


#set view 0, 0
# Plot the data
splot 'wavefunction' using 1:2:3 with points palette pointtype 7 pointsize 0.2
set output 'output2.png'


set view 90, 0

# Plot the data
splot 'wavefunction' using 1:2:3 with points palette pointtype 7 pointsize 0.2
set output 'output3.png'


set view 0, 0

# Plot the data
splot 'wavefunction' using 1:2:3 with points palette pointtype 7 pointsize 0.2
