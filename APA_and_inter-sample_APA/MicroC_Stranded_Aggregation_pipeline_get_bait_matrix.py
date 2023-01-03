import sys
import pandas as pd

'''
1 - contact (pairs) file output from distiler filtered using: zcat perfix.pairs.gz | awk 'BEGIN {OFS = "\t"} ; {if ($1 == "." & $2 == $4 & $9 >= 30 & $10 >= 30) {print $2, $3, $4, $5, $6, $7, $8, $9, $10}}' > perfix.nodups_30_intra.pairs
2 - region - the bait center
3 - baitStrand - the strand field for the bait
4 - bed2 file of regions of interest that you want to use as prey (the script will center the heatmap over the center of the regions)
5 - window half size of region to screen for contacts per region with other region in the list (integer)
6 - half the number of pixels to devide the window to (integer). The resulting pixels will represent a window-size/number-of-pixels bps
7 - matrix output (csv)
'''

contacts = pd.read_csv(sys.argv[1], names = ['chrA', 'sideA', 'chrB', 'sideB', 'strandA', 'strandB', 'read_type', 'mapqA', 'mapqB'], sep = '\t')
region = int(sys.argv[2])
baitStrand = str(sys.argv[3])
preys = pd.read_csv(sys.argv[4], names = ['chr', 'start', 'end', 'strand', 'center'], sep = '\t')
winSize = int(sys.argv[5])
pixNum = int(sys.argv[6])

pixelsList = list(range(-winSize, 0, int(winSize/pixNum))) + list(range(int(winSize/pixNum), winSize + 1, int(winSize/pixNum)))

# starting with an all zero dataframe
aggregate = pd.DataFrame(0, index = pixelsList, columns = pixelsList)
# This will result in each column representing a sub-region of the bait
# and each row representing a sub-region of the prey
for n in preys.index:
	regionB = preys['center'][n]
	if (baitStrand == '+') and (preys['strand'][n] == '+'):
		p = 0
		aggregatePair = pd.DataFrame(index = pixelsList, columns = pixelsList)
		a = region - winSize
		step = int(winSize/pixNum)
		while a < region + winSize:
			conRegA = contacts.loc[(contacts['sideA'].isin(list(range(a, a + step + 1)))) | (contacts['sideB'].isin(list(range(a, a + step + 1))))]
			addingList = []
			b = regionB - winSize
			while b < regionB + winSize: 
				conRegAB = conRegA.loc[(conRegA['sideA'].isin(list(range(b, b + step + 1))) & (conRegA['sideB'].isin(list(range(a, a + step + 1))))) | (conRegA['sideB'].isin(list(range(b, b + step + 1))) & (conRegA['sideA'].isin(list(range(a, a + step + 1)))))]
				addingList.append(len(conRegAB.index))
				b += step
			aggregatePair[pixelsList[p]] = addingList
			p += 1
			a += step
		aggregate = sum([aggregate, aggregatePair])
	elif (baitStrand == '+') and (preys['strand'][n] == '-'):
		p = 0
		aggregatePair = pd.DataFrame(index = pixelsList, columns = pixelsList)
		a = region - winSize
		step = int(winSize/pixNum)
		while a < region + winSize:
			conRegA = contacts.loc[(contacts['sideA'].isin(list(range(a, a + step + 1)))) | (contacts['sideB'].isin(list(range(a, a + step + 1))))]
			addingList = []
			b = regionB + winSize
			while b > regionB - winSize:
				conRegAB = conRegA.loc[(conRegA['sideA'].isin(list(range(b - step, b + 1))) & (conRegA['sideB'].isin(list(range(a, a + step + 1))))) | (conRegA['sideB'].isin(list(range(b - step, b + 1))) & (conRegA['sideA'].isin(list(range(a, a + step + 1)))))]
				addingList.append(len(conRegAB.index))
				b -= step
			aggregatePair[pixelsList[p]] = addingList
			p += 1
			a += step
		aggregate = sum([aggregate, aggregatePair])
	elif (baitStrand == '-') and (preys['strand'][n] == '+'):
		p = 0
		aggregatePair = pd.DataFrame(index = pixelsList, columns = pixelsList)
		a = region + winSize
		step = int(winSize/pixNum)
		while a > region - winSize:
			conRegA = contacts.loc[(contacts['sideA'].isin(list(range(a - step, a + 1)))) | (contacts['sideB'].isin(list(range(a - step, a + 1))))]
			addingList = []
			b = regionB - winSize
			while b < regionB + winSize:
				conRegAB = conRegA.loc[(conRegA['sideA'].isin(list(range(b, b + step + 1))) & (conRegA['sideB'].isin(list(range(a - step, a + 1))))) | (conRegA['sideB'].isin(list(range(b, b + step + 1))) & (conRegA['sideA'].isin(list(range(a - step, a + 1)))))]
				addingList.append(len(conRegAB.index))
				b += step
			aggregatePair[pixelsList[p]] = addingList
			p += 1
			a -= step
		aggregate = sum([aggregate, aggregatePair])
	elif (baitStrand == '-') and (preys['strand'][n] == '-'):
		p = 0
		aggregatePair = pd.DataFrame(index = pixelsList, columns = pixelsList)
		a = region + winSize
		step = int(winSize/pixNum)
		while a > region - winSize:
			conRegA = contacts.loc[(contacts['sideA'].isin(list(range(a - step, a + 1)))) | (contacts['sideB'].isin(list(range(a - step, a + 1))))]
			addingList = []
			b = regionB + winSize
			while b > regionB - winSize:
				conRegAB = conRegA.loc[(conRegA['sideA'].isin(list(range(b - step, b + 1))) & (conRegA['sideB'].isin(list(range(a - step, a + 1))))) | (conRegA['sideB'].isin(list(range(b - step, b + 1))) & (conRegA['sideA'].isin(list(range(a - step, a + 1)))))]
				addingList.append(len(conRegAB.index))
				b -= step
			aggregatePair[pixelsList[p]] = addingList
			p += 1
			a -= step
		aggregate = sum([aggregate, aggregatePair])

aggregate.iloc[::-1].to_csv(sys.argv[7])