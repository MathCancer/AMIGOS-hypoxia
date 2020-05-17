import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rk4 import rk4

def KI67_Basic ( t, rf, Par ):
  Q = rf[0]
  K = rf[1]
  
  r01 = Par[0] #1.0/(4.59*60.0)
  r10 = Par[1] #1.0/(15.5*60.0)
  
  dQdt =  -r01 * Q + 2*r10 * K
  dKdt = r01 * Q - r10 * K
  
  drfdt = np.array ( [ dQdt, dKdt] )

  return drfdt
  
def Model(Par):   
  #Po2 = np.array([4.60234, 4.74557, 4.95134, 5.26767, 5.70256, 6.22612, 6.87626, 7.68088, 8.67263, 9.85707, 11.2488, 12.8964, 14.8848, 17.2725, 20.1098, 23.3787, 27.324, 32.1182, 37.6101, 45.94])
  Po2 = np.array([0.0864051, 0.0890906, 0.0929483, 0.0988798, 0.107035, 0.116854, 0.129048, 0.144141, 0.162745, 0.184965, 0.211073, 0.241984, 0.279285, 0.324073, 0.377286, 0.438583, 0.512567, 0.602552, 0.706005, 0.831655, 0.97609, 1.15126, 1.3635, 1.60932, 1.90144, 2.25161, 2.67429, 3.17616, 3.76633, 4.46569, 5.30917, 6.32467, 7.53771, 8.94839, 10.6579, 12.746, 15.1701, 18.1281, 21.5589, 25.7389, 30.8287, 36.7742, 45.94])
  tspan = np.array ( [ 0.0, 6000.0 ] )
  InitialCond = np.array ( [ 56827, 0.0] )
  n = 300
  QOI1 = np.zeros((len(Po2),n+1))
  QOI2 = np.zeros((len(Po2),n+1))
  ParLocal = np.copy(Par)
  for index in range( 0,len(Po2) ):
    rate = (Po2[index] - ParLocal[3])/(ParLocal[2] - ParLocal[3])
    if(rate < 0.0): rate = 0.0
    if(rate > 1.0): rate = 1.0
    ParLocal[0] = Par[0]*rate
    t, ODE_out = rk4 ( KI67_Basic, tspan, InitialCond, n, ParLocal )
    QOI1[index,:] = 100*ODE_out[:,0]/(ODE_out[:,0]+ODE_out[:,1])
    QOI2[index,:] = 100*ODE_out[:,1]/(ODE_out[:,0]+ODE_out[:,1])
  QOI2 = np.delete(QOI2,-5,0)
  return QOI2[:,-1]

def ModelComplete(Par):   
  #Po2 = np.array([4.60234, 4.74557, 4.95134, 5.26767, 5.70256, 6.22612, 6.87626, 7.68088, 8.67263, 9.85707, 11.2488, 12.8964, 14.8848, 17.2725, 20.1098, 23.3787, 27.324, 32.1182, 37.6101, 45.94])
  Po2 = np.array([0.0864051, 0.0890906, 0.0929483, 0.0988798, 0.107035, 0.116854, 0.129048, 0.144141, 0.162745, 0.184965, 0.211073, 0.241984, 0.279285, 0.324073, 0.377286, 0.438583, 0.512567, 0.602552, 0.706005, 0.831655, 0.97609, 1.15126, 1.3635, 1.60932, 1.90144, 2.25161, 2.67429, 3.17616, 3.76633, 4.46569, 5.30917, 6.32467, 7.53771, 8.94839, 10.6579, 12.746, 15.1701, 18.1281, 21.5589, 25.7389, 30.8287, 36.7742, 45.94])
  tspan = np.array ( [ 0.0, 6000.0 ] )
  InitialCond = np.array ( [ 56827, 0.0] )
  n = 300
  QOI1 = np.zeros((len(Po2),n+1))
  QOI2 = np.zeros((len(Po2),n+1))
  ParLocal = np.copy(Par)
  for index in range( 0,len(Po2) ):
    rate = (Po2[index] - ParLocal[3])/(ParLocal[2] - ParLocal[3])
    if(rate < 0.0): rate = 0.0
    if(rate > 1.0): rate = 1.0
    ParLocal[0] = Par[0]*rate
    #print("Oxygen: "+str(Po2[index])+" Par: " + str(ParLocal[0])+" "+ str(Par[1])+" "+ str(Par[2])+" "+ str(Par[3]))
    t, ODE_out = rk4 ( KI67_Basic, tspan, InitialCond, n, ParLocal )
    #ODE_out[ODE_out > 1.0] = 1.0
    QOI1[index,:] = 100.0*ODE_out[:,0]/(ODE_out[:,0]+ODE_out[:,1])
    QOI2[index,:] = 100.0*ODE_out[:,1]/(ODE_out[:,0]+ODE_out[:,1])
  return QOI1,QOI2,t,Po2
  
#Data
with open('CalibProl.dat', 'r') as f:
  P1 = []
  P2 = []
  P3 = []
  P4 = []
  Ind = []
  AcT = []
  distance = []
  for each in f:
    each1 = each.split()
    P1.append(float(each1[0]))
    P2.append(float(each1[1]))
    P3.append(float(each1[2]))
    P4.append(float(each1[3]))
    Ind.append(int(each1[4]))
    AcT.append(int(each1[5]))
    distance.append(float(each1[6]))

Par1 = np.array(P1)
Par2 = np.array(P2)
Par3 = np.array(P3)
Par4 = np.array(P4)

#Plotting
# figure, axes = plt.subplots(nrows=1, ncols=3)
# plt.subplot(131)
# n, bins, patches = plt.hist(x=Par1, bins='auto', color='#009900',
              # alpha=0.7, rwidth=0.85)
# plt.grid(axis='y', alpha=0.75)
# plt.xlabel('$Persistence$ $time$')
# maxfreq1 = n.max()
# plt.ylim(ymax=np.ceil(maxfreq1 / 10) * 10 if maxfreq1 % 10 else maxfreq1 + 10)
# ProbPar1 = bins[np.argmax(n)]

# plt.subplot(132)
# n, bins, patches = plt.hist(x=Par2, bins='auto', color='#009900',
              # alpha=0.7, rwidth=0.85)
# plt.grid(axis='y', alpha=0.75)
# plt.xlabel('$Bias$')
# maxfreq2 = n.max()
# plt.ylim(ymax=np.ceil(maxfreq2 / 10) * 10 if maxfreq2 % 10 else maxfreq2 + 10)
# #plt.xlim((0.0, 0.035)) 
# ProbPar2 = bins[np.argmax(n)]

# plt.subplot(133)
# n, bins, patches = plt.hist(x=Par3, bins='auto', color='#009900',
              # alpha=0.7, rwidth=0.85)
# plt.grid(axis='y', alpha=0.75)
# plt.xlabel('$Speed$')
# maxfreq2 = n.max()
# plt.ylim(ymax=np.ceil(maxfreq2 / 10) * 10 if maxfreq2 % 10 else maxfreq2 + 10)
# #plt.xlim((0.0, 0.035)) 
# ProbPar3 = bins[np.argmax(n)]
# plt.suptitle("MAP = %2.4f ; %2.4f ; %2.4f" % (ProbPar1, ProbPar2, ProbPar3),fontsize=12)
# figure.tight_layout(pad=2.0)
# plt.savefig('HistPar.png')

figure, axes = plt.subplots(nrows=2, ncols=2,figsize=(12,12))

sns.set()
sns.set_color_codes("deep")
plt.subplot(221)
#sns.distplot(Par1)
value1 = sns.distplot(Par1,color='blue').get_lines()[0].get_data()
maxPar1 = value1[0][np.argmax(value1[1])]
plt.ylabel('Density')
plt.title("MAP = %2.4f and Min = %2.4f" % (maxPar1,Par1[np.argmin(distance)]), fontsize=12)
plt.xlabel('$r_{01}$')

plt.subplot(222)
sns.set_style('darkgrid')
#sns.distplot(Par2,color='blue',norm_hist=True)
value2 = sns.distplot(Par2,color='blue').get_lines()[0].get_data()
maxPar2 = value2[0][np.argmax(value2[1])]
plt.title("MAP = %2.4f and Min = %2.4f" % (maxPar2,Par2[np.argmin(distance)]), fontsize=12)
plt.xlabel('$r_{10}$')

plt.subplot(223)
sns.set_style('darkgrid')
#sns.distplot(Par3,color='blue',norm_hist=True)
value3 = sns.distplot(Par3,color='blue').get_lines()[0].get_data()
maxPar3 = value3[0][np.argmax(value3[1])]
plt.ylabel('Density')
plt.title("MAP = %2.4f and Min = %2.4f" % (maxPar3,Par3[np.argmin(distance)]), fontsize=12)
plt.xlabel('$\sigma_{S}$')

plt.subplot(224)
sns.set_style('darkgrid')
#sns.distplot(Par4,color='blue',norm_hist=True)
value4 = sns.distplot(Par4,color='blue').get_lines()[0].get_data()
maxPar4 = value4[0][np.argmax(value4[1])]
plt.title("MAP = %2.4f and Min = %2.4f" % (maxPar4,Par4[np.argmin(distance)]), fontsize=12)
plt.xlabel('$\sigma_{T}$')
figure.tight_layout(pad=1.0)
plt.savefig('Ki67.png')

plt.clf()
N = 20
fig = plt.figure()

DataCol1, DataCol2, DataCol3 = [], [], []
for line in open('../DataProl.dat', 'r'):
  values = [float(s) for s in line.split()]
  DataCol1.append(values[0])
  DataCol2.append(values[1])
  DataCol3.append(values[2])
DataCol1 = np.array(DataCol1)
DataCol2 = 55*np.array(DataCol2)
DataCol3 = 55*np.array(DataCol3)

#Par = np.array([maxPar1,maxPar2,maxPar3,maxPar4])
Par = np.array([0.01,0.01,16.5,0.0])


plt.errorbar(DataCol1[-N:], DataCol2[-N:], yerr=DataCol3[-N:],fmt='x')
plt.plot(DataCol1,Model(Par),color='blue')
plt.plot(DataCol1,Model([Par1[np.argmin(distance)],Par2[np.argmin(distance)],Par3[np.argmin(distance)],Par4[np.argmin(distance)]]),color='red')
plt.xlabel('Distance from the tumor core (mm)')
plt.ylabel('Ki67 (%)')
plt.savefig('Result.png')

#MAP Output
Out1, Out2, t, Po2 = ModelComplete(Par)
yr = np.linspace(t.min(), t.max(), len(t))
xr = np.linspace(0.1, 4.3, len(Po2))
X,Y = np.meshgrid(xr, yr)
figure, axes = plt.subplots(nrows=1, ncols=2,figsize=(12,6))
plt.subplot(121)
v = np.linspace(Out1.min(), Out1.max(), 10, endpoint=True)
plt.contourf(X,Y,Out1[:,:].T,v,cmap='binary');
x = plt.colorbar(ticks=v,label='Ki67- (%)')
plt.xlabel('Distance from the tumor core (mm)')
plt.ylabel('Time (min)')
plt.subplot(122)
v = np.linspace(Out2.min(), Out2.max(), 10, endpoint=True)
plt.contourf(X,Y,Out2[:,:].T,v,cmap='binary');
x = plt.colorbar(ticks=v,label='Ki67+ (%)')
plt.xlabel('Distance from the tumor core (mm)')
plt.ylabel('Time (min)')
figure.tight_layout(pad=1.0)
plt.savefig('Result2.png')  


plt.clf()
plt.xlabel('Distance from the tumor core (mm)')
plt.ylabel('$PO_2$ (mmHg)')
plt.plot(xr,Po2,color='blue')
plt.savefig('Result3.png')

#plt.figure(2)
#plt.ylim(Po2.min(), Po2.max())
#plt.xlim(t.min(), t.max())
#plt.title( 'Par1:'+str("%4.2f"%(Par[0]))+' - Par2:'+str("%4.2f"%(Par[1]))+ ' - Par3:' +str("%4.2f"%(Par[2])) + ' - Par4:' +str("%4.2f"%(Par[3])))

#row, col = Out1.shape
#for index in range( 0,row ):
# index = 40 
# plt.plot ( t, Out1[index,:], 'r-', linewidth = 2,color='r' )
# plt.plot ( t, Out2[index,:], 'r-', linewidth = 2,color='b' )
# plt.figure(2)

fig = plt.figure(4)
# plt.errorbar(X[-N:], Y[-N:], yerr=Ystd[-N:],ls='none')
# plt.title('Slice of the last time step')
plt.xlabel('Distance from the tumor core (mm)')
plt.ylabel('Ki67 (%)')
# plt.ylim(-0.01, 1.01)
plt.plot ( xr, Out1[:,-1], 'r-', linewidth = 2,color='r' )
plt.plot ( xr, Out2[:,-1], 'r-', linewidth = 2,color='b' )
plt.savefig('Result4.png')

#plt.plot ( xr, Out1[:,-50], '--', linewidth = 2,color='r' )
#plt.plot ( xr, Out2[:,-50], '--', linewidth = 2,color='b' )