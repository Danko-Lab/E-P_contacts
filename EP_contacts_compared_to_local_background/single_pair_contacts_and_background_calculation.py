import sys
import pandas as pd

'''
1 - promoter center
2 - enhancer center
3 - promoter contacts
4 - enhancer contacts
5 - winSize: half size of window arounf enhancers and promoters to count E-P contacts
6 - minBGdis: the minimum distance to take from element center for BG calculations
7 - maxBGdis: the maximum distance to take from element center for BG calculations
8 - output file
'''
promoter_center = int(sys.argv[1])
enhancer_center = int(sys.argv[2])

promoter_contacts = pd.read_csv(sys.argv[3], names = ['chrA', 'sideA', 'chrB', 'sideB', 'strandA', 'strandB', 'read_type', 'mapqA', 'mapqB'], sep = '\t')

enhancer_contacts = pd.read_csv(sys.argv[4], names = ['chrA', 'sideA', 'chrB', 'sideB', 'strandA', 'strandB', 'read_type', 'mapqA', 'mapqB'], sep = '\t')

winSize = int(sys.argv[5])
minBGdis = int(sys.argv[6])
maxBGdis = int(sys.argv[7])
# calculate enhancer-promoter contacts based on window size defined that sorrounds their centers
ep_contacts = promoter_contacts.loc[((promoter_contacts['sideA'].isin(list(range(promoter_center - winSize, promoter_center + winSize + 1)))) | (promoter_contacts['sideB'].isin(list(range(promoter_center - winSize, promoter_center + winSize + 1))))) & ((promoter_contacts['sideA'].isin(list(range(enhancer_center - winSize, enhancer_center + winSize + 1)))) | (promoter_contacts['sideB'].isin(list(range(enhancer_center - winSize, enhancer_center + winSize + 1)))))]

# calculate enhancer- and promoter-associated contacts near the other anchor, that are not enhancer-promoter contacts
promoter_contacts = promoter_contacts.loc[(promoter_contacts['sideA'].isin(list(range(enhancer_center - maxBGdis, enhancer_center - minBGdis + 1)) + list(range(enhancer_center + minBGdis, enhancer_center + maxBGdis + 1)))) | (promoter_contacts['sideB'].isin(list(range(enhancer_center - maxBGdis, enhancer_center - minBGdis + 1)) + list(range(enhancer_center + minBGdis, enhancer_center + maxBGdis + 1))))]
enhancer_contacts = enhancer_contacts.loc[(enhancer_contacts['sideA'].isin(list(range(promoter_center - maxBGdis, promoter_center - minBGdis + 1)) + list(range(promoter_center + minBGdis, promoter_center + maxBGdis + 1)))) | (enhancer_contacts['sideB'].isin(list(range(promoter_center - maxBGdis, promoter_center - minBGdis + 1)) + list(range(promoter_center + minBGdis, promoter_center + maxBGdis + 1))))]

#define EP contacts and background based on earlier calculations
ep = len(ep_contacts.index)
background = len(promoter_contacts.index) + len(enhancer_contacts.index)
#append the result for the pair to the results file
with open(sys.argv[8], "a") as file_object:
	file_object.write(promoter_contacts['chrA'][list(promoter_contacts.index)[0]] + '\t' + str(promoter_center) + '\t' + str(enhancer_center) + '\t' + str(ep) + '\t' + str(background) + '\n')