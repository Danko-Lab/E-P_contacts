### This script is an example of visualization of changes in
### enhancer-promoter contacts over local background ratio.
### The data used for this example is from the 30 minutes dTAG degredation
### of NELFB in mESCs compared to control (T=0).

from scipy.stats import gaussian_kde
import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#load the calculated enhancer-promoter and background contacts numbers 
#calculated and calculate ratios between them 
t0 = pd.read_csv(os.path.join('C:\\Users','barsh','Desktop','Gilad','Focal_Loops','mESCs_dTAG_T0_EP_and_BG_contacts.txt'), sep = '\t', names = ['chr', 'promoter', 'enhancer', 'EP_contacts', 'BG_contacts'])
t0['ratio'] = t0['EP_contacts'] / t0['BG_contacts']
t30min = pd.read_csv(os.path.join('C:\\Users','barsh','Desktop','Gilad','Focal_Loops','mESCs_dTAG_T30min_EP_and_BG_contacts.txt'), sep = '\t', names = ['chr', 'promoter', 'enhancer', 'EP_contacts', 'BG_contacts'])
t30min['ratio'] = t30min['EP_contacts'] / t30min['BG_contacts']

#generate a merged DF for plotting
merged_dTAG = t0.merge(t30min, how = 'inner', on = ['chr', 'promoter', 'enhancer'])
merged_dTAG.columns = ['chr', 'promoter', 'enhancer', 'EP_contacts_t0', 'BG_contacts_t0', 'ratio_t0',
               'EP_contacts_t30min', 'BG_contacts_t30min', 'ratio_t30min']
#define a minimum baseline number of contacts
merged_dTAG = merged_dTAG.loc[merged_dTAG['EP_contacts_t0'] >= 8]

#Plot a density scatterplot based on the data
values = np.vstack([merged_dTAG["ratio_t0"], merged_dTAG["ratio_t30min"]])
kernel = gaussian_kde(values)(values)
fig, ax = plt.subplots(figsize=(6, 6))
sns.scatterplot(
    data=merged_dTAG,
    x="ratio_t0",
    y="ratio_t30min",
    c=kernel,
    cmap="jet",
    ax=ax,
)
plt.xscale('log')
plt.yscale('log')
plt.xlim(plt.ylim())

xpoints = ypoints = plt.xlim()
plt.plot(xpoints, ypoints, linestyle='--', color='k', lw=3, scalex=False, scaley=False)
plt.ylabel('T = 30 min', size = 20)
plt.xlabel('T = 0', size = 20)
plt.title( 'E-P Contacts Over Background', size = 20)
plt.yticks(size = 16)
plt.xticks(size = 16)
merged_dTAG['ROR'] = merged_dTAG['ratio_t0'] / merged_dTAG['ratio_t30min']
up_1 = len([x for x in list(merged_dTAG['ROR']) if x > 1])
up_2= len([x for x in list(merged_dTAG['ROR']) if x < 1])

#print the number of enhancer-promoter pairs in total,
#and the number of pairs that lost (T0) or gained (T30min)
#contacts following the treatment
print('all: ' + str(len(merged_dTAG['ROR'])))
print('T0: ' + str(up_1))
print('T30min: ' + str(up_2))
