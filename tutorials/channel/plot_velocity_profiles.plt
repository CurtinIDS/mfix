###################################################################################################
# Define various constants
###################################################################################################
pi = acos(-1.0)
beta = 20.0                                               # sloping angle (Deg.)
Umean = 1.0                                               # Mean velocity magnitude (m/s)

alpha = 90.0 - beta
alpha = alpha /180.0 * pi
###################################################################################################
# Analytical solution - No-slip wall
###################################################################################################
u_mag_analytic(x) =  1.5 * Umean *(1.0 - x**2)           
u_analytic(x) =  u_mag_analytic(x) * cos(alpha)
v_analytic(x) =  u_mag_analytic(x) * sin(alpha)
###################################################################################################
# Analytical solution - Free-slip wall
###################################################################################################
#u_mag_analytic(x) =  Umean            
#u_analytic(x) =  u_mag_analytic(x) * cos(alpha)
#v_analytic(x) =  u_mag_analytic(x) * sin(alpha)
###################################################################################################
# Common setup
###################################################################################################
# Set output terminal - postscript allows more options like superscript, subscript etc
set term postscript enhanced  24 #color 24
set size 0.75,1.5
# No legend (key)
# set nokey
set key graph 0.8, graph 0.95
set key spacing 4
# set key box linestyle 1
# Grid is on
set grid
###################################################################################################
# U-velocity profile
###################################################################################################
# Title
set title "U-velocity profile ( {/Symbol b} = 20 ^o )"
# X and Y labels and range
set xlabel "x* / 2H [-]"
set ylabel "U/U_m [-]"
set xrange [-1.0 : 1.0]
set yrange [0.0 : 2.0]
# Output file name
set output "Ug_B20_FOUP.ps"
# Plotting command
plot "Ug_B20_FOUP.dat"      using     1:2     with points  pointtype 6  pointsize 2       title "Numerical" , \
      u_analytic(x)                           with lines   linetype 1   linewidth 4       title "Analytical"
# Output file name
# Usage image magic to rotate and dump as jpg
! convert -crop 0x0 -rotate 90 -quality 100 Ug_B20_FOUP.ps Ug_B20_FOUP.jpg
# Use display to view the output
! rm -f Ug_B20_FOUP.ps
! display Ug_B20_FOUP.jpg &
###################################################################################################
# V-velocity profile
###################################################################################################
# Title
set title "V-velocity profile ( {/Symbol b} = 20 ^o )"
# X and Y labels and range
set xlabel "x* / 2H [-]"
set ylabel "V/U_m [-]"
set xrange [-1.0 : 1.0]
set yrange [0.0 : 2.0]
# Output file name
set output "Vg_B20_FOUP.ps"
# Plotting command
plot "Vg_B20_FOUP.dat"     using    1:2         with points  pointtype 6  pointsize 2       title "Numerical" , \
      v_analytic(x)                             with lines   linetype 1   linewidth 4       title "Analytical"
# Usage image magic to rotate and dump as jpg
! convert -crop 0x0 -rotate 90 -quality 100 Vg_B20_FOUP.ps Vg_B20_FOUP.jpg
# Use display to view the output
! rm -f Vg_B20_FOUP.ps
! display Vg_B20_FOUP.jpg &

