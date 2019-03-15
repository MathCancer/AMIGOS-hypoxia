import glob
import os

for fname in glob.glob("output_%08d_vasculature.mat"):
  cmd = "python plot_vasculature.py " + fname
  print(cmd)
  os.system(cmd)
