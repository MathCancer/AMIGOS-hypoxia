import numpy as np
import matplotlib.pyplot as plt
import shlex,subprocess
from subprocess import DEVNULL, STDOUT, check_call

def Model(Par):   
  #Po2 = np.array([4.60234, 4.74557, 4.95134, 5.26767, 5.70256, 6.22612, 6.87626, 7.68088, 8.67263, 9.85707, 11.2488, 12.8964, 14.8848, 17.2725, 20.1098, 23.3787, 27.324, 32.1182, 37.6101, 45.94])
  Po2 = np.array([0.0864051, 0.0890906, 0.0929483, 0.0988798, 0.107035, 0.116854, 0.129048, 0.144141, 0.162745, 0.184965, 0.211073, 0.241984, 0.279285, 0.324073, 0.377286, 0.438583, 0.512567, 0.602552, 0.706005, 0.831655, 0.97609, 1.15126, 1.3635, 1.60932, 1.90144, 2.25161, 2.67429, 3.17616, 3.76633, 4.46569, 5.30917, 6.32467, 7.53771, 8.94839, 10.6579, 12.746, 15.1701, 18.1281, 21.5589, 25.7389, 30.8287, 36.7742, 45.94])
  QOI = np.ones(43)
  QOI = (Po2 - Par[1])/(Par[0] - Par[1])
  QOI[QOI < 0] = 0.0
  QOI[QOI > 1] = 1.0
  return QOI
  

def ABC_MCMC(Npar,data,std, LowLimit, UpperLimit, FILE='CalibMCMC.dat', tol = 66, max_iterations=10**10, NumAccept = 200):
  file = open(FILE,"w") 
  count = 0
  theta_star = np.zeros( (1,Npar) )
  for j in range(0, Npar):
    theta_star[0,j] = np.random.uniform(LowLimit[j],UpperLimit[j])

  for i in range(0, max_iterations):
    output_model = Model([theta_star[0,0],theta_star[0,1]])
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
with open('DataProl.dat', 'r') as f:
  KI67 = []
  KI67Std = []
  for each in f:
    each1 = each.split()
    KI67.append(float(each1[0]))
    KI67Std.append(float(each1[1]))

KI67 = np.array(KI67)
KI67Std = np.array(KI67Std)

UpperLimit = np.array([300.0,10.0])
LowLimit = np.array([45.0,0.0])
qoi = KI67
std = KI67Std

ABC_MCMC(UpperLimit.shape[0],qoi, std, LowLimit, UpperLimit)