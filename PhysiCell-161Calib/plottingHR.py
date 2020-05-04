from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt

initial_index = 0;
last_index = 336;
Lcell_size = 5;
Dcell_size = 1;
Ecell_size = 5;
live_count = np.zeros( last_index+1 );
dead_count = np.zeros( last_index+1 );
endo_count = np.zeros( last_index+1 );
times = np.zeros( last_index+1 );
for n in range( initial_index,last_index+1 ):
  filename='output'+"%08i"%n+'.xml'
  filenameOut='Figs/output'+"%08i"%n+'.png'
  mcds=pyMCDS(filename,'output')
  times[n]= mcds.get_time()
  
  cx = mcds.data['discrete_cells']['position_x'];
  cy = mcds.data['discrete_cells']['position_y'];
  cycle = mcds.data['discrete_cells']['cycle_model']
  cycle = cycle.astype( int )
  cell_type = mcds.data['discrete_cells']['cell_type']
  cell_type = cell_type.astype(int)
  live = np.argwhere( (cycle < 100) & (cell_type==0) ).flatten()
  dead = np.argwhere( (cycle >= 100) & (cell_type==0) ).flatten()
  endo = np.argwhere( cell_type==1 ).flatten()
  live_count[n] = len(live)
  dead_count[n] = len(dead)
  endo_count[n] = len(endo)

  # proteinsX = mcds.data['discrete_cells']['proteins_x']
  # proteinsY = mcds.data['discrete_cells']['proteins_y']
  # rgb_var = np.vstack((0.5*proteinsX[live],0.5*proteinsY[live],0*proteinsX[live])).T
	
  plt.clf()
  fig = plt.figure()
  o2 = mcds.get_concentrations( 'oxygen' );
  X,Y = mcds.get_2D_mesh();
  plt.contourf(X,Y,o2[:,:,0],cmap='binary');
  plt.colorbar()
  #plt.scatter( cx[live],cy[live],c=rgb_var,s=Lcell_size);
  plt.axis('image')
  plt.xlim(-1500, 1500)
  plt.ylim(-1500, 1500)
  plt.title( '#LC:'+str("%04i"%(live_count[n]))+'  #DC:'+str("%04i"%(dead_count[n]))+'  #EC:'+str("%04i"%endo_count[n])+ '  Time:' +str("%8.2f"%(n/24)) + ' days', size=10)
  #plt.scatter( cx[dead],cy[dead],c='moccasin',s=Dcell_size );
  #plt.scatter( cx[endo],cy[endo],c='r',s=Ecell_size );
  fig.savefig(filenameOut)


#rgb_var = np.vstack((0.5*proteinsX[live],0.5*proteinsY[live],0*proteinsX[live])).T
#print(len(cx))
#for i in range(len(cx)):
#  if cell_type[i] != 0:
#    print(cell_type[i]);

# plt.figure(1) 
# plt.clf()
# plt.scatter(cx[live],cy[live],c=rgb_var, s=Lcell_size);
# plt.scatter( cx[dead],cy[dead],c='moccasin',s=Dcell_size );
# plt.title( 'Cells colored by oncoprotein value' , size=10)
# plt.axis( 'image' )
# plt.xlim(-1500, 1500)
# plt.ylim(-1500, 1500)
# plt.xlabel( 'x' , size=10 )
# plt.ylabel( 'y', size=10 )


#plt.clf()
#plt.figure(2)
#plt.xlim(-1500, 1500)
#plt.ylim(-1500, 1500) 
#o2 = mcds.get_concentrations( 'oxygen' );
#vegf = mcds.get_concentrations( 'VEGF' );
#ecm = mcds.get_concentrations( 'ECM' );
#X,Y = mcds.get_2D_mesh();
#plt.contourf(X,Y,o2[:,:,0],cmap='binary');
#plt.colorbar()
#plt.axis( 'image' )

# plt.figure(3)
# plt.contourf(X,Y,vegf[:,:,0],cmap='binary');
# plt.colorbar()
# plt.axis( 'image' )

# plt.figure(4)
# plt.contourf(X,Y,ecm[:,:,0],cmap='binary');
# plt.colorbar()
# plt.axis( 'image' )


#plt.clf()
#plt.figure(5) 
#mcds.get_substrate_names();
#o2 = mcds.get_concentrations( 'oxygen' );
#X,Y = mcds.get_2D_mesh();
#plt.contourf(X,Y,o2[:,:,0],cmap='binary');
#plt.colorbar()
#plt.scatter( cx,cy,c='r',s=Lcell_size);
#plt.axis('image')
#plt.xlim(-1500, 1500)
#plt.ylim(-1500, 1500)
#plt.title( 'Live cells colored at t=' +str(t/60) + ' hr', size=10)
# let's plot dead cells as black
#plt.scatter( cx[dead],cy[dead],c='moccasin',s=Dcell_size );
#plt.scatter( cx[endo],cy[endo],c='r',s=Dcell_size );
#plt.show()
