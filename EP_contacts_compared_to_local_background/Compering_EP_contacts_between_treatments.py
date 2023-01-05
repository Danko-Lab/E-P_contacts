from scipy.stats import gaussian_kde
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

'''
1 - MicroC_EP_and_BG_contacts.bsh for DMSO control
2 - MicroC_EP_and_BG_contacts.bsh for FLV
3 - MicroC_EP_and_BG_contacts.bsh for TRP
4 - Contact depth for DMSO control
5 - Contact depth for FLV
6 - Contact depth for TRP
7 - minimum EP contacts in CPB
8 - output scatterplot for FLV vs. DMSO
9 - output scatterplot for TRP vs. DMSO
10 - output scatterplot for FLV vs. TRP

'''

dmso_contact_number = int(sys.argv[4])
flv_contact_number = int(sys.argv[5])
trp_contact_number = int(sys.argv[6])
ep_contact_threshold = int(sys.argv[7])

dmso = pd.read_csv(sys.argv[1], sep = '\t', names = ['chr', 'promoter', 'enhancer', 'EP_contacts', 'BG_contacts'])
dmso['ratio'] = dmso['EP_contacts'] / dmso['BG_contacts']
flv = pd.read_csv(sys.argv[2], sep = '\t', names = ['chr', 'promoter', 'enhancer', 'EP_contacts', 'BG_contacts'])
flv['ratio'] = flv['EP_contacts'] / flv['BG_contacts']
trp = pd.read_csv(sys.argv[3], sep = '\t', names = ['chr', 'promoter', 'enhancer', 'EP_contacts', 'BG_contacts'])
trp['ratio'] = trp['EP_contacts'] / trp['BG_contacts']
a = dmso.merge(trp, how = 'inner', on = ['chr', 'promoter', 'enhancer'])
merged_txn = a.merge(flv, how = 'inner', on = ['chr', 'promoter', 'enhancer'])
merged_txn.columns = ['chr', 'promoter', 'enhancer', 'EP_contacts_dmso', 'BG_contacts_dmso', 'ratio_dmso',
               'EP_contacts_trp', 'BG_contacts_trp', 'ratio_trp',
               'EP_contacts_flv', 'BG_contacts_flv', 'ratio_flv']

dmso_contacts = dmso_contact_number/10000000000
flv_contacts = flv_contact_number/10000000000
trp_contacts = trp_contact_number/10000000000


merged_txn['EP_CPB_dmso'] = merged_txn['EP_contacts_dmso']/dmso_contacts
merged_txn['EP_CPB_flv'] = merged_txn['EP_contacts_flv']/flv_contacts
merged_txn['EP_CPB_trp'] = merged_txn['EP_contacts_trp']/trp_contacts


merged_txn = merged_txn.loc[((merged_txn['EP_CPB_dmso'] > ep_contact_threshold) |
                            (merged_txn['EP_CPB_flv'] > ep_contact_threshold) |
                            (merged_txn['EP_CPB_trp'] > ep_contact_threshold)) &
                            ((merged_txn['EP_CPB_dmso'] > 0) &
                            (merged_txn['EP_CPB_flv'] > 0) &
                            (merged_txn['EP_CPB_trp'] > 0))].dropna()



values = np.vstack([merged_txn["ratio_dmso"], merged_txn["ratio_flv"]])
print(len([x for x in list(merged_txn["ratio_flv"]/merged_txn["ratio_dmso"]) if x > 1]))
print(len([x for x in list(merged_txn["ratio_flv"]/merged_txn["ratio_dmso"]) if x < 1]))
kernel = gaussian_kde(values)(values)
fig, ax = plt.subplots(figsize=(8, 8))
sns.scatterplot(
    data=merged_txn,
    x="ratio_dmso",
    y="ratio_flv",
    c=kernel,
    cmap="jet",
    ax=ax,
)
plt.xscale('log')
plt.yscale('log')
plt.xlim(plt.ylim())
xpoints = ypoints = plt.xlim()
plt.plot(xpoints, ypoints, linestyle='--', color='k', lw=3, scalex=False, scaley=False)
plt.ylabel('FLV', size = 30)
plt.xlabel('DMSO', size = 30)
plt.yticks(size = 20)
plt.xticks(size = 20)
plt.savefig(sys.argv[8])

values = np.vstack([merged_txn["ratio_dmso"], merged_txn["ratio_trp"]])
print(len([x for x in list(merged_txn["ratio_trp"]/merged_txn["ratio_dmso"]) if x > 1]))
print(len([x for x in list(merged_txn["ratio_trp"]/merged_txn["ratio_dmso"]) if x < 1]))
kernel = gaussian_kde(values)(values)
fig, ax = plt.subplots(figsize=(8, 8))
sns.scatterplot(
    data=merged_txn,
    x="ratio_dmso",
    y="ratio_trp",
    c=kernel,
    cmap="jet",
    ax=ax,
)
plt.xscale('log')
plt.yscale('log')
plt.xlim(plt.ylim())
xpoints = ypoints = plt.xlim()
plt.plot(xpoints, ypoints, linestyle='--', color='k', lw=3, scalex=False, scaley=False)
plt.ylabel('TRP', size = 30)
plt.xlabel('DMSO', size = 30)
plt.yticks(size = 20)
plt.xticks(size = 20)
plt.savefig(sys.argv[9])



values = np.vstack([merged_txn["ratio_trp"], merged_txn["ratio_flv"]])
print(len([x for x in list(merged_txn["ratio_flv"]/merged_txn["ratio_trp"]) if x > 1]))
print(len([x for x in list(merged_txn["ratio_flv"]/merged_txn["ratio_trp"]) if x < 1]))
kernel = gaussian_kde(values)(values)
fig, ax = plt.subplots(figsize=(6, 6))
sns.scatterplot(
    data=merged_txn,
    x="ratio_trp",
    y="ratio_flv",
    c=kernel,
    cmap="jet",
    ax=ax,
)
plt.xscale('log')
plt.yscale('log')
plt.xlim(plt.ylim())
xpoints = ypoints = plt.xlim()
plt.plot(xpoints, ypoints, linestyle='--', color='k', lw=3, scalex=False, scaley=False)
plt.ylabel('FLV', size = 20)
plt.xlabel('TRP', size = 20)
plt.yticks(size = 16)
plt.xticks(size = 16)
plt.savefig(sys.argv[10])
