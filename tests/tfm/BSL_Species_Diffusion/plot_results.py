#!/usr/bin/python

import numpy as np
import pylab as plt

# load data file and ignore top three text lines
with open('solution_xg_profile.dat') as f:
  data_lines = f.readlines()[3:]

data = np.loadtxt(data_lines)

# plot numerical and exact solution in first subplot
# plot error in second subplot
plt.figure(1, figsize=(9.0,4.0))

plt.subplot(1,2,1)

plot1 = plt.plot(data[:,0], data[:,3], 'o', mfc='none', mec='b', mew='1')
plot2 = plt.plot(data[:,0], data[:,5], '-', c='b', lw=2)
plot3 = plt.plot(data[:,0], data[:,4], 'o', mfc='none', mec='r', mew='1')
plot4 = plt.plot(data[:,0], data[:,6], '-', c='r', lw=2)
#plt.xlabel('x (m)',labelpad=10)
plt.xlabel('x (m)')
plt.ylabel('Species mass fraction')
plt.xlim(0.0, 0.5)
plt.ylim(0.0, 1.0)
plt.legend([plot1, plot2, plot3, plot4], ('X(Air):MFIX','X(Air): Exact',
  'X($H_2O$): MFIX','X($H_2O$): Exact'), 'best', prop={'size':6})
plt.title('(a)',x=0.5,y=-0.3)

plt.subplot(1,2,2)
plot5 = plt.plot(data[:,0], data[:,7], '-b+')
plot6 = plt.plot(data[:,0], data[:,8], '-r+')
#plt.xlabel('x (m)',labelpad=10)
plt.xlabel('x (m)')
plt.ylabel('Error')
#plt.locator_params(axis='x',nbins=5)
plt.legend([plot5, plot6], ('X($H_2O$)','X(Air)'), 'best',prop={'size':6})

plt.title('(b)',x=0.5,y=-0.3)

plt.subplots_adjust(wspace=0.3)

plt.savefig('tfm05_01.png', bbox_inches='tight')
plt.close(1)
