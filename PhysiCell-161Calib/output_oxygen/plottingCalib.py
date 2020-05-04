import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#Data
with open('OxyCalibMCMC.dat', 'r') as f:
  P1 = []
  P2 = []
  P3 = []
  Ind = []
  AcT = []
  distance = []
  for each in f:
    each1 = each.split()
    P1.append(float(each1[0]))
    P2.append(float(each1[1]))
    P3.append(float(each1[2]))
    Ind.append(int(each1[3]))
    AcT.append(int(each1[4]))
    distance.append(float(each1[5]))

Par1 = np.array(P1)
Par2 = np.array(P2)
Par3 = np.array(P3)

print("Min distance: "+str(distance[np.argmin(distance)])+" - Par:"+str(Par1[np.argmin(distance)])+" "+str(Par2[np.argmin(distance)])+" "+str(Par3[np.argmin(distance)])+"\n")

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

#figure, axes = plt.subplots(nrows=1, ncols=3,figsize=(12,6))

sns.set()
sns.set_color_codes("deep")
# plt.subplot(131)
# sns.distplot(Par1,color='blue')
# value1 = sns.distplot(Par1).get_lines()[0].get_data()
# maxPar1 = value1[0][np.argmax(value1[1])]
# plt.title("MAP = %2.4f" % maxPar1, fontsize=12)
# plt.xlabel('Persistence time')

# plt.subplot(132)
# sns.set_style('darkgrid')
# sns.distplot(Par2,color='blue')
# value2 = sns.distplot(Par2).get_lines()[0].get_data()
# maxPar2 = value2[0][np.argmax(value2[1])]
# plt.title("MAP = %2.4f" % maxPar2, fontsize=12)
# plt.xlabel('Bias')

# plt.subplot(133)
sns.set_style('darkgrid')
sns.distplot(Par3,color='blue')
value3 = sns.distplot(Par3).get_lines()[0].get_data()
maxPar3 = value3[0][np.argmax(value3[1])]
plt.title("MAP = %2.4f and Min = %2.4f" % (maxPar3,Par3[np.argmin(distance)]), fontsize=12)
plt.xlabel('Cell uptake')
plt.ylabel('Frequency')
#figure.tight_layout(pad=1.0)
plt.savefig('Oxygen.png')
