import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
'''
1 - ContactCaller_microC_output_W_functional_nonfunctional_and_oter_pair_assignments.txt
2 - Minimum contacts
3 - Minimum distance
4 - Output Violinplot
'''

contacts_data = pd.read_csv(sys.argv[1], sep = '\t', index_col = 0)
contacts_data.index = range(len(contacts_data.index))


contacts_data = contacts_data.loc[contacts_data['observed'] >= int(sys.argv[2])]

contacts_data = contacts_data.loc[contacts_data['expected'] >= int(sys.argv[2])]

contacts_data = contacts_data.loc[contacts_data['distance'] >= int(sys.argv[3])]

contacts_data['Obs/Exp'] = np.log2(contacts_data['observed'] / contacts_data['expected'])


plt.figure(figsize=(6,8))
sns.violinplot(x = 'PosNeg', y = 'Obs/Exp', data = contacts_data, showfliers = False, palette = ['#ffa600', '#bc5090', 'gray'], order = ['positive', 'negative', 'other'], inner = 'quartile')
plt.yticks(size = 16)
plt.ylabel('Normalized contacts (log2)', size = 24)
plt.xlabel('Pair Type', size = 24)
plt.xticks([0, 1, 2], ['Functional', 'Nonfunctional', 'Other'], size = 16)
plt.savefig(sys.argv[4])
plt.show()
