#!/usr/bin/python

import numpy as np
import pylab as plt

# load data files and ignore top three text lines
with open('u_profile_Re100_S.dat') as f:
  data_lines = f.readlines()[3:]

data_u_100 = np.loadtxt(data_lines)

with open('v_profile_Re100_S.dat') as f:
  data_lines = f.readlines()[3:]

data_v_100 = np.loadtxt(data_lines)

with open('u_profile_Re400_S.dat') as f:
  data_lines = f.readlines()[3:]

data_u_400 = np.loadtxt(data_lines)

with open('v_profile_Re400_S.dat') as f:
  data_lines = f.readlines()[3:]

data_v_400 = np.loadtxt(data_lines)

with open('u_profile_Re1000_S.dat') as f:
  data_lines = f.readlines()[3:]

data_u_1000 = np.loadtxt(data_lines)

with open('v_profile_Re1000_S.dat') as f:
  data_lines = f.readlines()[3:]

data_v_1000 = np.loadtxt(data_lines)

with open('ghia_results_u.dat') as f:
  data_lines = f.readlines()[3:]

data_ug = np.loadtxt(data_lines)

with open('ghia_results_v.dat') as f:
  data_lines = f.readlines()[3:]

data_vg = np.loadtxt(data_lines)


# plot vertical centerline comparisons in 1st sub plot
# plot horizontal centerline comparisons in 2nd sub plot

plt.figure(1, figsize=(9.0,8.0)) # size of the full figure in inches
# Re=100 comparisons:

plt.subplot(2,2,1)

plot1 = plt.plot(data_u_100[:,1], data_u_100[:,3], '-b', linewidth=2)
plot2 = plt.plot(data_ug[:,1], data_ug[:,2], 'ro')
#plt.xlabel('x (m)',labelpad=10)
plt.xlabel('y (m)')
plt.ylabel('U_g (m/s)')
#plt.xlim(0.0, 1.0)
#plt.ylim(0.0, 1.0)
plt.legend([plot1, plot2], ('MFIX','Ref.'), 'upper left')
plt.title('(a) x=0.5, Re=100',x=0.5,y=-0.25)

plt.subplot(2,2,2)
plot3 = plt.plot(data_v_100[:,0], data_v_100[:,3], '-b', linewidth=2)
plot4 = plt.plot(data_vg[:,0], data_vg[:,2], 'ro')
#plt.xlabel('x (m)',labelpad=10)
plt.xlabel('x (m)')
plt.ylabel('V_g (m/s)')
#plt.xlim(0.0, 1.0)
plt.ylim(-0.5, 0.4)
plt.legend([plot1, plot2], ('MFIX','Ref.'), 'best')
#plt.locator_params(axis='x',nbins=5)
plt.title('(b) y=0.5, Re=100',x=0.5,y=-0.25)

#plt.subplots_adjust(wspace=0.3)

#plt.savefig('tfm03_01.png', bbox_inches='tight')
#plt.close(1)


# Re=400 comparisons:
#plt.figure(1, figsize=(9.0,4.0)) # size of the full figure in inches

plt.subplot(2,2,3)

plot1 = plt.plot(data_u_400[:,1], data_u_400[:,3], '-b', linewidth=2)
plot2 = plt.plot(data_ug[:,1], data_ug[:,3], 'ro')
#plt.xlabel('x (m)',labelpad=10)
plt.xlabel('y (m)')
plt.ylabel('U_g (m/s)')
#plt.xlim(0.0, 1.0)
#plt.ylim(0.0, 1.0)
plt.legend([plot1, plot2], ('MFIX','Ref.'), 'upper left')
plt.title('(c) x=0.5, Re=400',x=0.5,y=-0.25)

plt.subplot(2,2,4)
plot3 = plt.plot(data_v_400[:,0], data_v_400[:,3], '-b', linewidth=2)
plot4 = plt.plot(data_vg[:,0], data_vg[:,3], 'ro')
#plt.xlabel('x (m)',labelpad=10)
plt.xlabel('x (m)')
plt.ylabel('V_g (m/s)')
#plt.xlim(0.0, 1.0)
plt.ylim(-0.5, 0.4)
plt.legend([plot1, plot2], ('MFIX','Ref.'), 'best')
#plt.locator_params(axis='x',nbins=5)
plt.title('(d) y=0.5, Re=400',x=0.5,y=-0.25)

plt.subplots_adjust(wspace=0.3)
plt.subplots_adjust(hspace=0.4)

#plt.show()

plt.savefig('tfm03_01.png', bbox_inches='tight')
plt.close(1)


plt.figure(2, figsize=(9.0,4.0)) # size of the full figure in inches
# Re=1000 comparisons:

plt.subplot(1,2,1)

plot1 = plt.plot(data_u_1000[:,1], data_u_1000[:,3], '-b', linewidth=2)
plot2 = plt.plot(data_ug[:,1], data_ug[:,4], 'ro')
#plt.xlabel('x (m)',labelpad=10)
plt.xlabel('y (m)')
plt.ylabel('U_g (m/s)')
#plt.xlim(0.0, 1.0)
#plt.ylim(0.0, 1.0)
plt.legend([plot1, plot2], ('MFIX','Ref.'), 'upper left')
plt.title('(a) x=0.5, Re=1000',x=0.5,y=-0.25)

plt.subplot(1,2,2)
plot3 = plt.plot(data_v_1000[:,0], data_v_1000[:,3], '-b', linewidth=2)
plot4 = plt.plot(data_vg[:,0], data_vg[:,4], 'ro')
#plt.xlabel('x (m)',labelpad=10)
plt.xlabel('x (m)')
plt.ylabel('V_g (m/s)')
#plt.xlim(0.0, 1.0)
#plt.ylim(-0.5, 0.4)
plt.legend([plot1, plot2], ('MFIX','Ref.'), 'best')
#plt.locator_params(axis='x',nbins=5)
plt.title('(b) y=0.5, Re=1000',x=0.5,y=-0.25)

plt.subplots_adjust(wspace=0.3)
plt.subplots_adjust(hspace=0.4)

#plt.show()

plt.savefig('tfm03_02.png', bbox_inches='tight')
plt.close(2)
