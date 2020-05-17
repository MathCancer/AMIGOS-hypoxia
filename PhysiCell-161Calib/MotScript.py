import numpy as np
import matplotlib.pyplot as plt
import shlex,subprocess
from subprocess import DEVNULL, STDOUT, check_call
from MCMC import ABC_MCMC

def ModelRed(Par):
  run = "./motility2D.exe 2 20.35 "+str(Par[0])+"  0.2769"
  #print(run)
  args = shlex.split(run)
  p = subprocess.Popen(args,stdout=DEVNULL)
  if p.wait() != 0:
    print("There was an error")
  output = np.loadtxt("./output0002", dtype='f', delimiter='\t')
  speed = np.array(output[:,1])
  return speed
  
def ModelGreen(Par):
  run = "./motility2D.exe 2 40.86 "+str(Par[0])+"  0.3760"
  #print(run)
  args = shlex.split(run)
  p = subprocess.Popen(args,stdout=DEVNULL)
  if p.wait() != 0:
    print("There was an error")
  output = np.loadtxt("./output0002", dtype='f', delimiter='\t')
  speed = np.array(output[:,1])
  return speed

#Data
Temp, DsRedSpeed, DsRedSpeedSTD, GFPSpeed, GFPSpeedSTD = [], [], [], [], []
for line in open('DataMot.dat', 'r'):
  values = [float(s) for s in line.split()]
  Temp.append(values[0])
  DsRedSpeed.append(values[1])
  DsRedSpeedSTD.append(values[2])
  GFPSpeed.append(values[3])
  GFPSpeedSTD.append(values[4])
Temp = np.array(Temp)
DsRedSpeed = np.array(DsRedSpeed)
DsRedSpeedSTD = np.array(DsRedSpeedSTD)
GFPSpeed = np.array(GFPSpeed)
GFPSpeedSTD = np.array(GFPSpeedSTD)

UpperLimit = np.array([1.0])
LowLimit = np.array([0.0])
qoi = DsRedSpeed

ABC_MCMC(ModelRed, qoi,LowLimit, UpperLimit,'CalibMotR.dat',10.0,1)
ABC_MCMC(ModelRed, qoi,LowLimit, UpperLimit,'CalibMotG.dat',10.0,1)
#input('Press ENTER to exit')
