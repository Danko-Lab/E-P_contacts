### This script is for calculating the total contacts mapped to each of the defined windows
### aroung baits\preys for a single chromosome

import sys
import pandas as pd

'''
1 - chromosome
2 - directory for: whole chromosome contact (pairs) file output from distiler filtered using: zcat perfix.pairs.gz | awk 'BEGIN {OFS = '\t'} ; {if ($1 == '.' & $2 == $4 & $9 >= 30 & $10 >= 30) {print $2, $3, $4, $5, $6, $7, $8, $9, $10}}' > perfix.nodups_30_intra.pairs
3 - directory for: regions - the bait in chrom file
4 - baits or preys - whether the TREs are the baits or the preys in the original pile-up analysis
5 - window half size of region to screen for contacts per region with other region in the list (integer)
6 - half the number of pixels to devide the window to (integer). The resulting pixels will represent a window-size/number-of-pixels bps
7 - matrix output (csv)
'''
chr = str(sys.argv[1])
treType = str(sys.argv[4])
winSize = int(sys.argv[5])
pixNum = int(sys.argv[6])
pixelsList = list(range(-winSize, 0, int(winSize/pixNum))) + list(range(int(winSize/pixNum), winSize + 1, int(winSize/pixNum)))

aggregate = pd.DataFrame(0, index = pixelsList, columns = ['contacts'])

#Get the chromosome contact file. If contacts were not shifted, change "shifted" below to the appropriate suffix
contacts = pd.read_csv(str(sys.argv[2]) + chr + '.shifted', names = ['chrA', 'sideA', 'chrB', 'sideB', 'strandA', 'strandB', 'read_type', 'mapqA', 'mapqB'], sep = '\t')

#Get the bed file of the baits\preys
regions = pd.read_csv(sys.argv[3] + chr + '_' + treType + '.bed', names = ['chr', 'start', 'end', 'strand', 'tsn'], sep = '\t')

# starting with an all zero dataframe
# This will result in a table with each row representing a sub-region of the bait\prey
for i in regions.index:
	contactsi = []
	if (regions['strand'][i] == '+'):
		p = 0
		aggregatePair = pd.DataFrame(index = pixelsList, columns = ['contacts'])
		a = regions['tsn'][i] - winSize
		step = int(winSize/pixNum)
		while a < regions['tsn'][i] + winSize:
			conRegA = contacts.loc[(contacts['sideA'].isin(list(range(a, a + step + 1)))) | (contacts['sideB'].isin(list(range(a, a + step + 1))))]
			p += 1
			a += step
			contactsi.append(len(conRegA.index))
		aggregatePair['contacts'] = contactsi
	elif (regions['strand'][i] == '-'):
		p = 0
		aggregatePair = pd.DataFrame(index = pixelsList, columns = ['contacts'])
		a = regions['tsn'][i] + winSize
		step = int(winSize/pixNum)
		while a > regions['tsn'][i] - winSize:
			conRegA = contacts.loc[(contacts['sideA'].isin(list(range(a - step, a + 1)))) | (contacts['sideB'].isin(list(range(a - step, a + 1))))]
			p += 1
			a -= step
			contactsi.append(len(conRegA.index))
		aggregatePair['contacts'] = contactsi
	aggregate = sum([aggregate, aggregatePair])
	print (chr + ' ' + treType + ' ' + str(i+1) + ' out of ' + str(len(regions.index)))
#Save per-chromosome global contacts table
aggregate.to_csv(sys.argv[7])