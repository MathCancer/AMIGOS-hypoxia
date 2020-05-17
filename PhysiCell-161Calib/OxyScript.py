import numpy as np
import matplotlib.pyplot as plt
import shlex,subprocess
from subprocess import DEVNULL, STDOUT, check_call
from MCMC import ABC_MCMC

def Model(Par):
  run = "./oxygen2D 1 0.1 100000.0 "+str(Par[0])
  #print(run)
  args = shlex.split(run)
  p = subprocess.Popen(args,stdout=DEVNULL)
  if p.wait() != 0:
    print("There was an error")
  output = np.loadtxt("./output0001", dtype='f', delimiter='\t') 
  oxygen = np.array(output[:,1])
  return oxygen

#Data
Dist, Oxygen, OxygenStd = [], [], []
for line in open('DataOxy.dat', 'r'):
  values = [float(s) for s in line.split()]
  Dist.append(values[0])
  Oxygen.append(values[1])
  OxygenStd.append(values[2])
Dist = np.array(Dist)
Oxygen = np.array(Oxygen)
OxygenStd = np.array(OxygenStd)

UpperLimit = np.array([2.0])
LowLimit = np.array([0.0])
qoi = Oxygen

ABC_MCMC(Model, qoi, LowLimit, UpperLimit,'CalibOxy.dat',30.0,1)