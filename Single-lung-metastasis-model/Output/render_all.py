import glob
import os

for fname in glob.glob("output_*.mat"):
  cmd = "python john2x4_abs_suptitle.py " + fname
  print(cmd)
  os.system(cmd)
