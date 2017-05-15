#!/usr/bin/python

import numpy as np
import pylab as plt

# load data file and ignore top three text lines
with open('solution_x_velocity_profile.dat') as f:
  data_lines = f.readlines()[3:]

data = np.loadtxt(data_lines)

# plot numerical and exact solution in first subplot
# plot error in second subplot
plt.figure(1, figsize=(7.0,3.0))

plt.subplot(1,2,1)

plot1 = plt.plot(data[:,3], data[:,1], 'ro')
plot2 = plt.plot(data[:,4], data[:,1], '-b+', linewidth=2)
plt.xlabel('x-velocity (m/s)',labelpad=20)
plt.ylabel('y (m)')
plt.xlim(0.0, 20.0)
plt.ylim(0.0, 0.01)
plt.legend([plot1, plot2], ('MFIX','Exact'), 'best')
plt.title('(a)',x=0.5,y=-0.4)

plt.subplot(1,2,2)
plot3 = plt.plot(data[:,5], data[:,1], '-b+')
plt.xlabel('Error (m/s)',labelpad=20)
plt.ylabel('y (m)')
plt.locator_params(axis='x',nbins=5)
plt.title('(b)',x=0.5,y=-0.4)

plt.subplots_adjust(wspace=0.5)

plt.savefig('tfm01_01.png', bbox_inches='tight')
plt.close(1)
