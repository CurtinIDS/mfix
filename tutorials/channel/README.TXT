This tutorial solve the gas flow in a channel, rotated 20 degrees with respect
to the vertical direction. 
The file usr3.f will extract velocity profiles at locations specified 
in extract_velocity.inp.
Assuming gnuplot is installed, the velocity profiles can be plotted using 
the command: gnuplot plot_velocity_profiles.plt
MFIX must be compiled from this directory to activate usr3.f

Note: To switch to Free-slip walls, 
Replace  BC_TYPE(3) = 'CG_NSW' by BC_TYPE(3) = 'CG_FSW' around line 124 in fmfix.dat
Replace NSW by FSW at line 1 of extract_velocity.inp
Uncomment the analytical solution calculation in plot_velocity_profiles.plt (lines 19-21)

