### This is a script for visualizing the inter-sample comperison of APAs

import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu
import scipy.stats as stats
import sys

'''
1 - Control observed contacts APA
2 - Treatment observed contacts APA
3 - Control baits 1D signal vector
4 - Control preys 1D signal vector
5 - Treatment baits 1D signal vector
6 - Treatment preys 1D signal vector
7 - winSize (half the 1D size of the APA (10000 for a 20Kb X 20Kb APA)
8 - pixNum (Half the number of cells to have in each row and each column of the APA)
9 - Number of baits (promoters in our case)
10 - Number of preys (enhancers or all TREs in our case)
11 - pathe to save the 1D normalized change APA matrics
'''



#Define the window size around anchors and pixel size for the APA
winSize = int(sys.argv[7])
pixNum = int(sys.argv[8])

#Define the window prey and bait numbers
baitNum = int(sys.argv[9])
preyNum = int(sys.argv[10])

# Load the APAs matrices calculatef for the different samples
control = pd.read_csv(sys.argv[1], index_col = 0)
control.columns = control.columns.astype(int)

treatment = pd.read_csv(sys.argv[2], index_col = 0)
treatment.columns = treatment.columns.astype(int)

# Load the tables of total contacts mapped to windows around baits and preys and devide by the number of baits\preys
#baits:
norcontrol = pd.read_csv(sys.argv[3], index_col = 0)/baitNum
nortreatment = pd.read_csv(sys.argv[5], index_col = 0)/baitNum

#preys:
norcontrolT = pd.read_csv(sys.argv[4], index_col = 0)/preyNum
nortreatmentT = pd.read_csv(sys.argv[6], index_col = 0)/preyNum

# calculate expected change matrices and obs/exp change matrices
lolET = []
lolT = []
for i in control.index:
    rowiET = []
    rowiT = []
    for j in control.columns:
        rowiET.append(np.sum([nortreatmentT['contacts'][i],nortreatment['contacts'][j]]) / np.sum([norcontrolT['contacts'][i],norcontrol['contacts'][j]]))
        rowiT.append(treatment[j][i]/float(control[j][i]))

    lolET.append(rowiET)
    lolT.append(rowiT)

expMatrixT = pd.DataFrame(lolET, index = control.index, columns = control.columns)
obsMatrixT = pd.DataFrame(lolT, index = control.index, columns = control.columns)
obsMatrixT.index = obsMatrixT.index.astype(int)
obsMatrixT.columns = obsMatrixT.columns.astype(int)
obsOVexpMatrixT = obsMatrixT.astype(float) / expMatrixT.astype(float)



#Heatmap presentation for treatment/control APA
ax = sns.heatmap(np.log2(obsOVexpMatrixT), cmap = "RdYlBu_r", square = True, vmax = 1, vmin = -1)
plt.xticks([0,pixNum, pixNum*2], ['-' + str(int(winSize / 1000)) + 'kb', '0', str(int(winSize / 1000)) + 'kb'])
plt.xlabel('Distane to Promoter TSS', size =16)
plt.yticks([0,pixNum, pixNum*2], [str(int(winSize / 1000)) + 'kb', '0', '-' + str(int(winSize / 1000)) + 'kb'])
plt.ylabel('Distane to Enhancer TSS', size =16)

plt.savefig(sys.argv[11])
