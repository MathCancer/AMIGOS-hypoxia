#
# Display scalar fields on a hexagonal mesh, using MATLAB formatted input files: 
# mesh.mat and output_*.mat
#
import matplotlib
matplotlib.use('Agg')
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scipy
from scipy import io
import numpy as np

print(len(sys.argv))
if (len(sys.argv) < 2):
  print("Usage: " + sys.argv[0] + " output_suffix.mat")
  sys.exit(0)
else:
  out_mat_file = sys.argv[1]

fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(14,6))
ax0, ax1, ax2, ax3, ax4, ax5, ax6, ax7 = axes.flatten()

print("----- mesh.mat ----")
d={}
scipy.io.loadmat("mesh",d)
m=d["mesh"]
print(m.shape)
x=m[1,:]
y=m[2,:]

print("X range=",x.min(),x.max())
print("Y range=",y.min(),y.max())

"""
    fwrite( (char*) &( voxel_population_vectors[i].live_cell_counts.at(0) ) , sizeof(double) , 1 , fp );// Live tumor cells
    fwrite( (char*) &( voxel_population_vectors[i].apoptotic_cell_counts.at(0) ) , sizeof(double) , 1 , fp ); //  apoptotic tumor cells
    fwrite( (char*) &( voxel_population_vectors[i].necrotic_cell_counts[0] ) , sizeof(double) , 1 , fp );  //  necrotic tumor cells
    fwrite( (char*) &( voxel_population_vectors[i].live_cell_counts.at(1) ) , sizeof(double) , 1 , fp );  //  Vascular cells
    fwrite( (char*) &( voxel_population_vectors[i].apoptotic_cell_counts.at(1) ) , sizeof(double) , 1 , fp );  //  Vasculars - NOT BEING USED - don't visualize me unless you just need to.  I would rather have a workign script that produces a blank animation/chart over nothing.
    fwrite( (char*) &( voxel_population_vectors[i].necrotic_cell_counts[1] ) , sizeof(double) , 1 , fp );   //  Vasular field - NOT BEING USED
    fwrite( (char*) &( substrate_vectors[i].substrate_quantity.at(0) ) , sizeof(double) , 1 , fp );  //  Oxygen
    fwrite( (char*) &( substrate_vectors[i].substrate_quantity.at(1) ) , sizeof(double) , 1 , fp );  //  Anigogenic Factor
    fwrite( (char*) &(  t ) , sizeof(double) , 1 , fp );  // time in minutes
"""
print("----- voxels (output*.mat) ----")
d2 = {}
scipy.io.loadmat(out_mat_file,d2)
v = d2["voxels_populations"]
print(v.shape)

# Use this block to give varying scales

#print("range of Live=",v[0,:].min(),v[0,:].max())
#print("range of Vascular=",v[3,:].min(),v[3,:].max())
#print("range of v[5,]=",v[5,:].min(),v[5,:].max())
#print("range of oxygen=",v[6,:].min(),v[6,:].max())
#print("range of angiogenic=",v[7,:].min(),v[7,:].max())

# Use this block to give constant color scales

print("range of Live=",0,29) # Set for hexoganal lattic with voxel diamter of 40 um
print("range of Vascular=",0,1)
print("range of v[5,]=",v[5,:].min(),v[5,:].max())
print("range of oxygen=",0,1)
print("range of angiogenic=",0,1)

live=v[0,:]
vasc=v[3,:]
o2=v[6,:]
angio=v[7,:]
time=v[8,1]
#print (time)
#-----------------------
# Options for custom colormaps: 
# https://matplotlib.org/examples/pylab_examples/custom_cmap.html
colors = [(0, 0, 1), (0, 1, 0), (1, 1, 0), (1, 0, 0)] 
cmap_name = 'my_colors'
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=4)  # N= # of discrete colors 

#------ Oxygen ---------
marker_size = 120  # area of (hexagonal) glyphs used for all scatter plots

half_dist = 300
ax0.axis([-half_dist,half_dist, -half_dist,half_dist])

# https://matplotlib.org/api/markers_api.html
# https://matplotlib.org/devdocs/api/_as_gen/matplotlib.pyplot.scatter.html
#matplotlib.pyplot.scatter(x, y, s=None, c=None, marker=None, cmap=None, norm=None, vmin=None, vmax=None, alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None, **kwargs)
voxel_plot0 = ax0.scatter(x, y, marker='h', s=marker_size, c=o2, cmap=cm, vmin=0.0, vmax=1.0)
ax0.set_aspect('equal')
ax0.set_title('Oxygen')
fig.colorbar(voxel_plot0, ax=ax0)

#------- Angiogenic factor --------
ax1.axis([-half_dist,half_dist, -half_dist,half_dist])

angio_min = 0
angio_max = 1
angio_del = (angio_max - angio_min)/3.0

voxel_plot1 = ax1.scatter(x, y, marker='h', s=marker_size, c=angio, cmap=cm, vmin=0.0, vmax=1.0)
ax1.set_aspect('equal')
ax1.set_title('Angio')
fig.colorbar(voxel_plot1, ax=ax1)

#------- Vasculature --------
ax2.axis([-half_dist,half_dist, -half_dist,half_dist])
voxel_plot2 = ax2.scatter(x, y, marker='h', s=marker_size, c=vasc, cmap=cm, vmin=0.0, vmax=1.0)
ax2.set_aspect('equal')
ax2.set_title('Vascular')
fig.colorbar(voxel_plot2, ax=ax2)

#------- Live cells --------
ax3.axis([-half_dist,half_dist, -half_dist,half_dist])
voxel_plot3 = ax3.scatter(x, y, marker='h', s=marker_size, c=live, cmap=cm, vmin=0.0, vmax=29.3)
ax3.set_aspect('equal')
ax3.set_title('Live')
fig.colorbar(voxel_plot3, ax=ax3)

#=================================================
#------- Histograms --------

# https://matplotlib.org/devdocs/api/_as_gen/matplotlib.pyplot.hist.html
#matplotlib.pyplot.hist(x, bins=None, range=None, density=None, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=None, log=False, color=None, label=None, stacked=False, normed=None, hold=None, data=None, **kwargs)
ax4.set_title('Histogram of O2')
n, bins, patches = ax4.hist(o2, 10, range=(0.0,1.0))

ax5.set_title('Histogram of Angio')
n, bins, patches = ax5.hist(angio, 10, range=(0.0,1.0))

ax6.set_title('Histogram of Vascular')
n, bins, patches = ax6.hist(v[3,:], 10, range=(0.0,1.0))

ax7.set_title('Histogram of Live')
n, bins, patches = ax7.hist(live, 10, range=(0.0,30.0))

#-----------------
plt.suptitle(time, y=0.999) # closer - needs to be more visible - may need to change canvas size or title location.  Then need to get the time actually out ...
fig.tight_layout()
# https://matplotlib.org/devdocs/api/_as_gen/matplotlib.pyplot.suptitle.html
#matplotlib.pyplot.suptitle(*args, **kwargs)

#plt.show()
png_file = out_mat_file[:-4]+".png"
print("---> " + png_file)
fig.savefig(png_file)
#fig.savefig(png_file, transparent=True)
