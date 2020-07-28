import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append('../')
from pyMCDS import pyMCDS

def ReadDataR():
  file = np.loadtxt("DispRed1.dat", dtype='f', delimiter='\t') 
  col1 = np.array(file[:,0])
  col2 = np.array(file[:,1])
  new1 = np.reshape(col1, (61, -1),order='F') 
  new2 = np.reshape(col2, (61, -1),order='F') 
  return new1,new2
  
def ReadDataG():
  file = np.loadtxt("DispGreen1.dat", dtype='f', delimiter='\t') 
  col1 = np.array(file[:,0])
  col2 = np.array(file[:,1])
  new1 = np.reshape(col1, (61, -1),order='F') 
  new2 = np.reshape(col2, (61, -1),order='F') 
  return new1,new2

def ReadSimulationR():
  initial_index = 0;
  last_index = 60;
  positionx = np.empty((0, 100), float);
  positiony = np.empty((0, 100), float);
  
  for n in range( initial_index,last_index+1 ):
    filename='output'+"%08i"%n+'.xml'
    mcds=pyMCDS(filename,'../outputMotRed')
  
    cx = mcds.data['discrete_cells']['position_x'];
    cy = mcds.data['discrete_cells']['position_y'];
    positionx = np.append(positionx,np.array([cx]),axis=0)
    positiony = np.append(positiony,np.array([cy]),axis=0)
  return positionx,positiony
  
def ReadSimulationG():
  initial_index = 0;
  last_index = 60;
  positionx = np.empty((0, 100), float);
  positiony = np.empty((0, 100), float);
  
  for n in range( initial_index,last_index+1 ):
    filename='output'+"%08i"%n+'.xml'
    mcds=pyMCDS(filename,'../outputMotGre')
  
    cx = mcds.data['discrete_cells']['position_x'];
    cy = mcds.data['discrete_cells']['position_y'];
    positionx = np.append(positionx,np.array([cx]),axis=0)
    positiony = np.append(positiony,np.array([cy]),axis=0)
  return positionx,positiony

def PlotResultR(posx, posy,posxDat,posyDat):
  theta = np.linspace(0, 2*np.pi, 100)
  r1 = 50
  r2 = 100
  r3 = 150
  r4 = 200
  r5 = 250
  #sns.boxplot(data=0);
  figure, axes = plt.subplots(nrows=1, ncols=2,figsize=(12,6))
  plt.subplot(121)
  plt.grid(linestyle='dashdot',linewidth=0.7,alpha = 0.3)
  plt.plot(r1*np.cos(theta), r1*np.sin(theta),'--',color='gray',alpha=1.0,linewidth=0.7)
  plt.plot(r2*np.cos(theta), r2*np.sin(theta),'--',color='gray',alpha=1.0,linewidth=0.7)
  plt.plot(r3*np.cos(theta), r3*np.sin(theta),'--',color='gray',alpha=1.0,linewidth=0.7)
  plt.plot(r4*np.cos(theta), r4*np.sin(theta),'--',color='gray',alpha=1.0,linewidth=0.7)
  plt.plot(r5*np.cos(theta), r5*np.sin(theta),'--',color='gray',alpha=1.0,linewidth=0.7)
  print("Number of trajectories (data): "+str(posxDat.shape))
  posx = np.delete(posx,np.s_[posxDat.shape[1]::],1)
  posy = np.delete(posy,np.s_[posyDat.shape[1]::],1)
  print("Number of trajectories (model): "+str(posx.shape))
  for i in range( 0,posxDat.shape[1]-1 ):
    plt.plot(posxDat[:,i],posyDat[:,i],color='gray',ls='solid',linewidth=2.0, alpha=1.0)
  plt.plot(posxDat[:,-1],posyDat[:,-1],color='gray',ls='solid',linewidth=2.0, alpha=1.0,label='Data')
  for i in range( 0,posx.shape[1]-1 ):
    plt.plot(posx[:,i],posy[:,i],color='red',ls='solid',linewidth=2.0, alpha=0.5)
  plt.plot(posx[:,-1],posy[:,-1],color='red',ls='solid',linewidth=2.0, alpha=0.5,label='Model')
  plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='center', ncol=2,)
  plt.text(270, -300, '$\mu m$', fontsize=10) 
  plt.text(-320, 270, '$\mu m$', fontsize=10)   
  #plt.axes().set_aspect(1)
  plt.xlim(-280.0,280.0)
  plt.xticks(np.arange(-250, 255, step=100.0))
  plt.yticks(np.arange(-250, 255, step=100.0))
  
  sns.set()
  sns.set_style("white")
  DistData = np.sqrt(posxDat[-1,:]**2 +  posyDat[-1,:]**2)
  DistModel =  np.sqrt(posx[-1,:]**2 +  posy[-1,:]**2) 
  plt.subplot(122)
  sns.distplot(DistData,color='Gray',bins=range(0, int(max(DistData))+25, 25),kde=False)
  sns.distplot(DistModel,color='red',bins=range(0, int(max(DistData))+25, 25),kde=False)
  plt.xlabel('Displacement')
  plt.ylabel('Frequency')
  
  plt.show()

def PlotResultG(posx, posy,posxDat,posyDat):
  theta = np.linspace(0, 2*np.pi, 100)
  r1 = 50
  r2 = 100
  r3 = 150
  r4 = 200
  r5 = 250
  #sns.boxplot(data=0);
  figure, axes = plt.subplots(nrows=1, ncols=2,figsize=(12,6))
  plt.subplot(121)
  plt.grid(linestyle='dashdot',linewidth=0.7,alpha = 0.3)
  plt.plot(r1*np.cos(theta), r1*np.sin(theta),'--',color='gray',alpha=1.0,linewidth=0.7)
  plt.plot(r2*np.cos(theta), r2*np.sin(theta),'--',color='gray',alpha=1.0,linewidth=0.7)
  plt.plot(r3*np.cos(theta), r3*np.sin(theta),'--',color='gray',alpha=1.0,linewidth=0.7)
  plt.plot(r4*np.cos(theta), r4*np.sin(theta),'--',color='gray',alpha=1.0,linewidth=0.7)
  plt.plot(r5*np.cos(theta), r5*np.sin(theta),'--',color='gray',alpha=1.0,linewidth=0.7)
  print("Number of trajectories (data): "+str(posxDat.shape))
  posx = np.delete(posx,np.s_[posxDat.shape[1]::],1)
  posy = np.delete(posy,np.s_[posyDat.shape[1]::],1)
  print("Number of trajectories (model): "+str(posx.shape))
  for i in range( 0,posxDat.shape[1]-1 ):
    plt.plot(posxDat[:,i],posyDat[:,i],color='gray',ls='solid',linewidth=2.0, alpha=1.0)
  plt.plot(posxDat[:,-1],posyDat[:,-1],color='gray',ls='solid',linewidth=2.0, alpha=1.0,label='Data')
  for i in range( 0,posx.shape[1]-1 ):
    plt.plot(posx[:,i],posy[:,i],color='green',ls='solid',linewidth=2.0, alpha=0.5)
  plt.plot(posx[:,-1],posy[:,-1],color='green',ls='solid',linewidth=2.0, alpha=0.5,label='Model')
  plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='center', ncol=2,)
  plt.text(270, -300, '$\mu m$', fontsize=10) 
  plt.text(-320, 270, '$\mu m$', fontsize=10)   
  #plt.axes().set_aspect(1)
  plt.xlim(-280.0,280.0)
  plt.xticks(np.arange(-250, 255, step=100.0))
  plt.yticks(np.arange(-250, 255, step=100.0))
  
  sns.set()
  sns.set_style("white")
  DistData = np.sqrt(posxDat[-1,:]**2 +  posyDat[-1,:]**2)
  DistModel =  np.sqrt(posx[-1,:]**2 +  posy[-1,:]**2) 
  plt.subplot(122)
  sns.distplot(DistData,color='Gray',bins=range(0, int(max(DistData))+25, 25),kde=False)
  sns.distplot(DistModel,color='green',bins=range(0, int(max(DistData))+25, 25),kde=False)
  plt.xlabel('Displacement')
  plt.ylabel('Frequency')
  
  plt.show()
  

# Red cells
posx, posy = ReadSimulationR()
posxDat,posyDat = ReadDataR()
PlotResultR(posx, posy ,posxDat,posyDat)

# Green cells
posx, posy = ReadSimulationG()
posxDat,posyDat = ReadDataG()
PlotResultG(posx, posy ,posxDat,posyDat)