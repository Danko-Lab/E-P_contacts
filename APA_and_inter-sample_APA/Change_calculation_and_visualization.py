### This is a script for visualizing the inter-sample comperison of APAs
### This example was used to compare FLV- and TRP-treate mESC samples 
### to DMSO-treated control. Data was taken from Hsieh et al. 2020 (GSE130275)
### DMSO - GSM3735106
### TRP - GSM3735107, GSM4173519, GSM4173520
### FLV - GSM3735108, GSM4173521, GSM4173522

import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu
import scipy.stats as stats



#Define the window size around anchors and pixel size for the APA
winSize = 10000
pixNum = 50

# Load the APAs matrices calculatef for the different samples
dmso = pd.read_csv(os.path.join('C:\\Users','barsh','Desktop','Gilad','Focal_Loops','AggMat_10kb_pixNum50_mESCs_DMSO.csv'), index_col = 0)
dmso.columns = dmso.columns.astype(int)

trp = pd.read_csv(os.path.join('C:\\Users','barsh','Desktop','Gilad','Focal_Loops','AggMat_10kb_pixNum50_mESCs_TRP.csv'), index_col = 0)
trp.columns = trp.columns.astype(int)

flv = pd.read_csv(os.path.join('C:\\Users','barsh','Desktop','Gilad','Focal_Loops','AggMat_10kb_pixNum50_mESCs_FLV.csv'), index_col = 0)
flv.columns = flv.columns.astype(int)

# Load the tables of total contacts mapped to windows around baits and preys and devide by the number of baits\preys
#baits:
#In this example, 10,599 baits were defined in the dataset
norDMSO = pd.read_csv(os.path.join('C:\\Users','barsh','Desktop','Gilad','Focal_Loops','All_bites_within_search_window_normalization_mESCs_DMSO.csv'), index_col = 0)/10599.0
norTRP = pd.read_csv(os.path.join('C:\\Users','barsh','Desktop','Gilad','Focal_Loops','All_bites_within_search_window_normalization_mESCs_TRP.csv'), index_col = 0)/10599.0
norFLV = pd.read_csv(os.path.join('C:\\Users','barsh','Desktop','Gilad','Focal_Loops','All_bites_within_search_window_normalization_mESCs_FLV.csv'), index_col = 0)/10599.0
#preys:
#In this example, 45,677 preys were defined in the dataset

norDMSOT = pd.read_csv(os.path.join('C:\\Users','barsh','Desktop','Gilad','Focal_Loops','All_preys_within_search_window_normalization_mESCs_DMSO.csv'), index_col = 0)/45677.0
norTRPT = pd.read_csv(os.path.join('C:\\Users','barsh','Desktop','Gilad','Focal_Loops','All_preys_within_search_window_normalization_mESCs_TRP.csv'), index_col = 0)/45677.0
norFLVT = pd.read_csv(os.path.join('C:\\Users','barsh','Desktop','Gilad','Focal_Loops','All_preys_within_search_window_normalization_mESCs_FLV.csv'), index_col = 0)/45677.0

# calculate expected change matrices and obs/exp change matrices
lolET = []
lolT = []
lolEF = []
lolF = []
for i in dmso.index:
    rowiET = []
    rowiT = []
    rowiEF = []
    rowiF = []
    for j in dmso.columns:
        rowiET.append((norTRPT['contacts'][i] + norTRP['contacts'][j]) / (norDMSOT['contacts'][i] + norDMSO['contacts'][j]))
        rowiT.append(trp[j][i]/float(dmso[j][i]))
        rowiEF.append((norFLVT['contacts'][i] + norFLV['contacts'][j]) / (norDMSOT['contacts'][i] + norDMSO['contacts'][j]))
        rowiF.append(flv[j][i]/float(dmso[j][i]))

    lolET.append(rowiET)
    lolT.append(rowiT)
    lolEF.append(rowiEF)
    lolF.append(rowiF)

expMatrixT = pd.DataFrame(lolET, index = dmso.index, columns = dmso.columns)
obspMatrixT = pd.DataFrame(lolT, index = dmso.index, columns = dmso.columns)
obspMatrixT.index = obspMatrixT.index.astype(int)
obspMatrixT.columns = obspMatrixT.columns.astype(int)
obsOVexppMatrixT = obspMatrixT.astype(float) / expMatrixT.astype(float)

expMatrixF = pd.DataFrame(lolEF, index = dmso.index, columns = dmso.columns)
obspMatrixF = pd.DataFrame(lolF, index = dmso.index, columns = dmso.columns)
obspMatrixF.index = obspMatrixF.index.astype(int)
obspMatrixF.columns = obspMatrixF.columns.astype(int)
obsOVexppMatrixF = obspMatrixF.astype(float) / expMatrixF.astype(float)


#Heatmap presentation for TRP/DMSO APA
ax = sns.heatmap(np.log2(obsOVexppMatrixT), cmap = "RdYlBu_r", square = True, vmax = 1, vmin = -1)
plt.xticks([0,pixNum, pixNum*2], ['-' + str(int(winSize / 1000)) + 'kb', '0', str(int(winSize / 1000)) + 'kb'])
plt.yticks([0,pixNum, pixNum*2], [str(int(winSize / 1000)) + 'kb', '0', '-' + str(int(winSize / 1000)) + 'kb'])
plt.show()
#Smoothed heatmap presentation for TRP/DMSO APA
lowerRes = {}
c = 1
index = []
while c < len(obsOVexppMatrixT.columns):
    lowerRes[obsOVexppMatrixT.columns[c]] = []    
    i = 1
    while i < len(obsOVexppMatrixT.index):
        lowerRes[obsOVexppMatrixT.columns[c]].append(np.mean(obsOVexppMatrixT[list(obsOVexppMatrixT.columns)[c-1:c+2]].loc[list(obsOVexppMatrixT.index)[i-1:i+2]].mean()))
        i += 1
    c += 1
obsOVexppMatrixTLR = pd.DataFrame(lowerRes, index = obsOVexppMatrixT.index[1::1]) 

sns.heatmap(np.log2(obsOVexppMatrixTLR), cmap = "RdYlBu_r", square = True, vmin = -1,vmax = 1)
plt.xticks([0,49.5, 99], ['-' + str(int(winSize / 1000)) + 'kb', 'TSN', str(int(winSize / 1000)) + 'kb'], size = 18, rotation = 0)
plt.yticks([0,49.5, 99], [str(int(winSize / 1000)) + 'kb', 'TSN', '-' + str(int(winSize / 1000)) + 'kb'], size = 18)
plt.show()

#Heatmap presentation for FLV/DMSO APA
ax = sns.heatmap(np.log2(obsOVexppMatrixF), cmap = "RdYlBu_r", square = True, vmax = 1, vmin = -1)
plt.xticks([0,pixNum, pixNum*2], ['-' + str(int(winSize / 1000)) + 'kb', '0', str(int(winSize / 1000)) + 'kb'])
plt.yticks([0,pixNum, pixNum*2], [str(int(winSize / 1000)) + 'kb', '0', '-' + str(int(winSize / 1000)) + 'kb'])
plt.show()

#Smoothed heatmap presentation for FLV/DMSO APA
lowerRes = {}
c = 1
index = []
while c < len(obsOVexppMatrixF.columns):
    lowerRes[obsOVexppMatrixF.columns[c]] = []    
    i = 1
    while i < len(obsOVexppMatrixF.index):
        lowerRes[obsOVexppMatrixF.columns[c]].append(np.mean(obsOVexppMatrixF[list(obsOVexppMatrixF.columns)[c-1:c+2]].loc[list(obsOVexppMatrixF.index)[i-1:i+2]].mean()))
        i += 1
    c += 1
obsOVexppMatrixFLR = pd.DataFrame(lowerRes, index = obsOVexppMatrixF.index[1::1]) 

sns.heatmap(np.log2(obsOVexppMatrixFLR), cmap = "RdYlBu_r", square = True, vmin = -1,vmax = 1)
plt.xticks([0,49.5, 99], ['-' + str(int(winSize / 1000)) + 'kb', 'TSN', str(int(winSize / 1000)) + 'kb'], size = 18, rotation = 0)
plt.yticks([0,49.5, 99], [str(int(winSize / 1000)) + 'kb', 'TSN', '-' + str(int(winSize / 1000)) + 'kb'], size = 18)
plt.show()

#Plotting box-plots for signal distribution in APA's dots, stripes and edges
#dots are defined as a 2kbX2kb window at the center of the APA matrix.
#stripe for each direction is defined as 4kbX8kb window, startin 2kb from the center of the APA matrix.
#edges are all 4kbX4kb windows at the corners of the APA matrix.
dot_vs_stipe_GB = pd.DataFrame(columns = ['mean_contacts', 'CID', 'CIU', 'W', 'P', 'treatments', 'Position'])
matrixList = [obsOVexppMatrixF, obsOVexppMatrixT]
treatmentList = ['FLV', 'TRP']
for i in range(len(matrixList)):
    matrixList[i].columns = matrixList[i].columns.astype(int)
    dot = matrixList[i].loc[[x for x in list(matrixList[i].index) if -1000<=x<=1000]][[x for x in list(matrixList[i].columns) if -1000<=x<=1000]]
    stripeL = matrixList[i].loc[[x for x in list(matrixList[i].index) if -2000<=x<=2000]][[x for x in list(matrixList[i].columns) if ((-2000>x))]]
    stripeR = matrixList[i].loc[[x for x in list(matrixList[i].index) if -2000<=x<=2000]][[x for x in list(matrixList[i].columns) if ((x>2000))]]
    stripeU = matrixList[i].loc[[x for x in list(matrixList[i].index) if ((x>2000))]][[x for x in list(matrixList[i].columns) if -2000<=x<=2000]].transpose()
    stripeD = matrixList[i].loc[[x for x in list(matrixList[i].index) if ((-2000>x))]][[x for x in list(matrixList[i].columns) if -2000<=x<=2000]].transpose()
    edges = pd.concat([matrixList[i].loc[[x for x in list(matrixList[i].index) if ((-6000>x) or (x>6000))]][[x for x in list(matrixList[i].columns) if ((-6000>x) or (x>6000))]], matrixList[i].loc[[x for x in list(matrixList[i].index) if ((-6000>x) or (x>6000))]][[x for x in list(matrixList[i].columns) if ((-6000>x) or (x>6000))]].transpose()])
   
    fig = plt.figure(figsize =(10, 7)) 
    ax = fig.add_subplot(111) 

    bp = ax.boxplot([dot.stack(), stripeU.stack(), stripeD.stack(), stripeR.stack(), stripeL.stack(), edges.stack()],
                showmeans = True, showfliers = False, labels = ['Dot', 'Upper\nStripe', 'Lower\nStripe', 'Right\nStripe', 'Left\nStripe', 'Edges'],
                meanprops={"marker":"s","markerfacecolor":"white", "markeredgecolor":"black"},
                patch_artist=True, notch = True)
    for patch in bp['boxes']:
        patch.set_facecolor('gray') 
        patch.set_linewidth(3) 
    
    for median in bp['medians']: 
        median.set(color ='red', linewidth = 3) 
    
    for whisker in bp['whiskers']:
        whisker.set(color ='black', linewidth = 3) 
    for cap in bp['caps']: cap.set(color ='black', linewidth = 3) 
    
    
    plt.yticks(size = 24)
    plt.xticks(size = 24)
    
    plt.ylim(top = 3)
    plt.ylabel('Observe/Expected Change', size = 24)
    s12 = 'ns'
    if mannwhitneyu(dot.stack(), stripeU.stack())[1] < 0.05:
        s12 = '*'
        if mannwhitneyu(dot.stack(), stripeU.stack())[1] < 0.001:
            s12 = '**'
        if mannwhitneyu(dot.stack(), stripeU.stack())[1] < 0.0001:
            s12 = '***'
    
    x1, x2 = 1, 2
    y, col = 1.4, 'k'
    plt.plot([x1, x1, x2, x2], [y, y, y, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y, s12, ha='center', va='bottom', color=col, size = 20)
    
    s13 = 'ns'
    if mannwhitneyu(dot.stack(), stripeD.stack())[1] < 0.05:
        s13 = '*'
        if mannwhitneyu(dot.stack(), stripeD.stack())[1] < 0.001:
            s13 = '**'
        if mannwhitneyu(dot.stack(), stripeD.stack())[1] < 0.0001:
            s13 = '***'
    
    x1, x2 = 1, 3
    y, col = 1.53, 'k'
    plt.plot([x1, x1, x2, x2], [y, y, y, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y, s13, ha='center', va='bottom', color=col, size = 20)
    
    s14 = 'ns'
    if mannwhitneyu(dot.stack(), stripeR.stack())[1] < 0.05:
        s14 = '*'
        if mannwhitneyu(dot.stack(), stripeR.stack())[1] < 0.001:
            s14 = '**'
        if mannwhitneyu(dot.stack(), stripeR.stack())[1] < 0.0001:
            s14 = '***'
    
    x1, x2 = 1, 4
    y, col = 1.63, 'k'
    plt.plot([x1, x1, x2, x2], [y, y, y, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y, s14, ha='center', va='bottom', color=col, size = 20)
    
    s15 = 'ns'
    if mannwhitneyu(dot.stack(), stripeL.stack())[1] < 0.05:
        s15 = '*'
        if mannwhitneyu(dot.stack(), stripeL.stack())[1] < 0.001:
            s15 = '**'
        if mannwhitneyu(dot.stack(), stripeL.stack())[1] < 0.0001:
            s15 = '***'
    
    x1, x2 = 1, 5
    y, col = 1.75, 'k'
    plt.plot([x1, x1, x2, x2], [y, y, y, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y, s15, ha='center', va='bottom', color=col, size = 20)
    
    s23 = 'ns'
    if mannwhitneyu(stripeU.stack(), stripeD.stack())[1] < 0.05:
        s23 = '*'
        if mannwhitneyu(stripeU.stack(), stripeD.stack())[1] < 0.001:
            s23 = '**'
        if mannwhitneyu(stripeU.stack(), stripeD.stack())[1] < 0.0001:
            s23 = '***'
    
    x1, x2 = 2, 3
    y, col = 1.45, 'k'
    plt.plot([x1, x1, x2, x2], [y, y, y, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y, s23, ha='center', va='bottom', color=col, size = 20)
    
    s45 = 'ns'
    if mannwhitneyu(stripeR.stack(), stripeL.stack())[1] < 0.05:
        s45 = '*'
        if mannwhitneyu(stripeR.stack(), stripeL.stack())[1] < 0.001:
            s45 = '**'
        if mannwhitneyu(stripeR.stack(), stripeL.stack())[1] < 0.0001:
            s45 = '***'
    
    x1, x2 = 4, 5
    y, col = 1.65, 'k'
    plt.plot([x1, x1, x2, x2], [y, y, y, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y, s45, ha='center', va='bottom', color=col, size = 20)
    
    s56 = 'ns'
    if mannwhitneyu(stripeL.stack(), edges.stack())[1] < 0.05:
        s56 = '*'
        if mannwhitneyu(stripeL.stack(), edges.stack())[1] < 0.001:
            s56 = '**'
        if mannwhitneyu(stripeL.stack(), edges.stack())[1] < 0.0001:
            s56 = '***'
    
    x1, x2 = 5, 6
    y, col = 1.69, 'k'
    plt.plot([x1, x1, x2, x2], [y, y, y, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+0.01, s56, ha='center', va='bottom', color=col, size = 20)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.ylim(0.2, 1.8)
    plt.yticks([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8])
    plt.show()

