import numpy as np
from random import random

def Generate2DInitialCond():
    tumor_radius = 2000
    cell_radius = 8.413
    volume_total = 20.0*np.pi*tumor_radius*tumor_radius
    volume_experiment = 50.0*347.0*347.0
    numCell_experim = 674
    factor = volume_experiment/numCell_experim
    numCell_sample = round(volume_total/factor)
    numCell_temp = 0

    rad = tumor_radius*random()
    theta = 2*np.pi*random()
    position = np.array([[rad*np.cos(theta),rad*np.sin(theta),0.0]])

    print(numCell_sample)

    while numCell_temp < numCell_sample:
      include = True
      rad = tumor_radius*random()
      theta = 2*np.pi*random()
      sample = np.array([rad*np.cos(theta),rad*np.sin(theta),0.0])
      for i in range(len(position)):
        d = np.linalg.norm(position[i,:] - sample, axis=0)
        if (d < 0.95*2*cell_radius):
          include = False
          break
      if (include == True):
        position = np.append(position,np.reshape(sample, (1, 3)),axis=0)
        numCell_temp += 1
      
    with open('Posfile.txt','w') as f:
      np.savetxt(f, position, fmt='%e')

def Generate3DInitialCond():
    tumor_radius = 400
    cell_radius = 8.413
    volume_total = (4.0/3.0)*np.pi*tumor_radius*tumor_radius*tumor_radius
    volume_experiment = 50.0*347.0*347.0
    numCell_experim = 674
    density = numCell_experim/volume_experiment
    numCell_sample = round(volume_total*density)
    numCell_temp = 0

    rad = tumor_radius*random()
    theta = np.pi*random()
    phi = 2*np.pi*random()
    position = np.array([[rad*np.sin(theta)*np.cos(phi),rad*np.sin(theta)*np.sin(phi),rad*np.cos(theta)]])

    print(numCell_sample)

    while numCell_temp < numCell_sample:
      include = True
      rad = tumor_radius*random()
      theta = np.pi*random()
      phi = 2*np.pi*random()
      sample = np.array([rad*np.sin(theta)*np.cos(phi),rad*np.sin(theta)*np.sin(phi),rad*np.cos(theta)])
      for i in range(len(position)):
        d = np.linalg.norm(position[i,:] - sample, axis=0)
        if (d < 0.95*2*cell_radius):
          include = False
          break
      if (include == True):
        position = np.append(position,np.reshape(sample, (1, 3)),axis=0)
        numCell_temp += 1
      
    with open('Posfile3D.txt','w') as f:
      np.savetxt(f, position, fmt='%e')

Generate3DInitialCond()