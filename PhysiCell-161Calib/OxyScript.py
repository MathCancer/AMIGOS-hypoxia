import numpy as np
import matplotlib.pyplot as plt
import shlex,subprocess
from subprocess import DEVNULL, STDOUT, check_call

def Model(Par):
  run = "./oxygen2D 1 "+str(Par[0])+" "+str(Par[1])+" "+str(Par[2])
  #print(run)
  args = shlex.split(run)
  p = subprocess.Popen(args,stdout=DEVNULL)
  if p.wait() != 0:
    print("There was an error")

  output = np.loadtxt("./output0001", dtype='f', delimiter='\t') 
  oxygen = np.array(output[:,1])
  return oxygen
  

def ABC_MCMC(Npar,data,std, LowLimit, UpperLimit, FILE='CalibMCMC.dat', tol = 1.0, max_iterations=10**10, NumAccept = 200):
  file = open(FILE,"w") 
  count = 0
  theta_star = np.zeros( (1,Npar) )
  for j in range(0, Npar):
    theta_star[0,j] = np.random.uniform(LowLimit[j],UpperLimit[j])

  for i in range(0, max_iterations):
    output_model = Model([theta_star[0,0],theta_star[0,1],theta_star[0,2]])
    distance = np.sqrt(np.sum([(a - b)**2/(c*c) for a, b, c in zip(output_model, data, std)]))
    if (distance < tol or count == 0):
        theta_star1 = theta_star
        count = count + 1
        for j in range(0, Npar):
          file.write(str(theta_star[0,j])+" ")
        file.write(str(count)+" "+str(i)+" "+str(distance)+"\n")
        print(str(count)+"/"+str(i)+" -- distance: "+str(distance)+"\n")

    if (count == NumAccept):
        break
    cond = True
    #print(distance)
    while(cond):
      noise = noise = np.random.normal(0, 0.2*(UpperLimit-LowLimit))
      theta_star = theta_star1 + noise
      cond = [False for k in range(0,Npar) if theta_star[0,k]>UpperLimit[k] or theta_star[0,k]<LowLimit[k]]
  file.close()

#Data
with open('DataOxy.dat', 'r') as f:
  Oxygen = []
  OxygenStd = []
  for each in f:
    each1 = each.split()
    Oxygen.append(float(each1[0]))
    OxygenStd.append(float(each1[1]))

Oxygen = np.array(Oxygen)
OxygenStd = np.array(OxygenStd)

UpperLimit = np.array([0.1,100000.0,2.0])
LowLimit = np.array([0.1,100000.0,0.0])
qoi = Oxygen
std = OxygenStd

ABC_MCMC(UpperLimit.shape[0],qoi, std, LowLimit, UpperLimit)