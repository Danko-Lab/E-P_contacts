import pandas as pd
import numpy as np
import sys

'''
1 - ContactCaller_microC_output
2 - CSV file with functional EP pairs based on CRISPRi
3 - CSV file with functional EP pairs based on CRISPRi
4 - path to ContactCaller_microC_output with functiona\nonfunctiona\other assignment
'''


# get contacts data and add distance
contacts_data = pd.read_csv(sys.argv[1], sep = '\t', usecols = [0,1,2,3,5,6], names = ['target_site.chr', 'target_site.center', 'target_promoter.center', 'directional_distance', 'observed', 'expected'])
contacts_data['distance'] = abs(contacts_data['directional_distance'])

# get the files of positive and negative sets of pairs from Gasperini
positives = pd.read_csv(sys.argv[2])
negatives = pd.read_csv(sys.argv[3])

# Assign pairs to positive, negative and other
posindex = []
negindex = []

for i in contacts_data.index:
    pairList = [contacts_data['target_site.chr'][i], contacts_data['target_promoter.center'][i], contacts_data['target_site.center'][i]]
    # get indices of positives and negatives
    for p in positives.index:
        if pairList == [positives['target_site.chr'][p], positives['target_promoter.center'][p], positives['target_site.center'][p]]:
            posindex.append(i)
            print('pos')
    for n in negatives.index:
        if pairList == [negatives['target_site.chr'][n], negatives['target_promoter.center'][n], negatives['target_site.center'][n]]:
            negindex.append(i)
            print('neg')

    print ('index_assignment', i, ' out of ', len(contacts_data.index))

posneg_column = []
for i in contacts_data.index:
    if i in posindex:
        posneg_column.append('positive')
    elif i in negindex:
        posneg_column.append('negative')
    else:
        posneg_column.append('other')
    print ('getting posNeg column', i, ' out of ', len(contacts_data.index))

contacts_data['PosNeg'] = posneg_column



contacts_data.to_csv(sys.argv[4], sep = '\t')
