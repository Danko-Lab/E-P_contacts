import sys
import pandas as pd
'''
1 - the baits\preys bed file
2 - directory for: files of contacts around promoters max TSNs
3 - baits or preys - whether the TREs are the baits or the preys in the original pile-up analysis
4 - window half size of region to screen for contacts per region with other region in the list (integer)
5 - half the number of pixels to devide the window to (integer). The resulting pixels will represent a window-size/number-of-pixels bps
6 - matrix output (csv)
'''
treType = str(sys.argv[3])
winSize = int(sys.argv[4])
pixNum = int(sys.argv[5])
pixelsList = list(range(-winSize, 0, int(winSize/pixNum))) + list(range(int(winSize/pixNum), winSize + 1, int(winSize/pixNum)))

aggregate = pd.DataFrame(0, index = pixelsList, columns = ['contacts'])

regions = pd.read_csv(sys.argv[1], names = ['chr', 'start', 'end', 'strand'], sep = '\t')

regions['tsn'] = ((regions['start'] + regions['end'])/2).astype(int)

# starting with an all zero dataframe
# This will result in each each row representing a sub-region of the bait
for i in regions.index:
	contacts = pd.read_csv(sys.argv[2] + regions['chr'][i] + '_' + str(regions['tsn'][i]) + '.bed', names = ['chrA', 'sideA', 'chrB', 'sideB', 'strandA', 'strandB', 'type', 'QA', 'QB'], sep = '\t')
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
	print (treType + ' ' + str(i+1) + ' out of ' + str(len(regions.index)))

#Save the vector of 1D signal associated with bait\prey set of loci
aggregate.to_csv(sys.argv[6])