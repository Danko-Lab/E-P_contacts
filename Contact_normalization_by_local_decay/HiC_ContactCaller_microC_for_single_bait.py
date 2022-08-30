# Note - This script also tests for a the significance of changes in observed
# Vs. expected changes in contacts using Fisher's exact test

import sys
import pandas as pd
import numpy as np
import statsmodels.api as sm
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import IntVector as ivect
import rpy2.robjects as robjects
import rpy2
import gzip

mass = importr('MASS')
stats = importr('stats')
base = importr('base')

peakFile = sys.argv[1] # File of bait peaks to analyze (format: Chromosome <\t> bait-center-position)
preyFile =  sys.argv[2] # specific bait's prey file to analyze (format: Chromosome <\t> prey-center-position)
conFile = sys.argv[3] # Contact file in juicer short format (given from merged_nodup.txt)
out = sys.argv[4] # Path to output file
DIST = int(sys.argv[5]) # Half size of "prey" search window around "bait" position
CAP = int(sys.argv[6]) # Half size of "bait" and "prey" contact capture window
minDist = 5000 #int(sys.argv[7]) # Minimum distance between bait and prey

# Calculate observed and expected distributions for each locus
plus, minus, contactProbabilities = [],[],[]
locus = preyFile.split('/')[1].split('.')[0].split('_')
chrom, position = locus[0], int(locus[1])

peakFile_DF = pd.read_csv(peakFile, names = ['chrom', 'pos'], sep = '\t')
peakFile_DF = peakFile_DF.loc[(peakFile_DF['chrom'] == chrom) & (peakFile_DF['pos'] == position)]
chrom = chrom.strip("chr")

preyFile_DF = pd.read_csv(preyFile, names = ['chrom', 'pos'], sep = '\t')

# read the contacts around loci from the splited temp_contact_file
models = int(DIST/50000)
bgModel = []


with gzip.open(conFile+"_IN_"+peakFile+"_folder/temp_"+str(list(peakFile_DF.index)[0])+".gz") as contact_fp:
	contacts= np.array([l.strip().split() for l in contact_fp.readlines()])
	print ("contacts with locus number " + str(list(peakFile_DF.index)[0]))
	contacts = contacts[:,1:4]
	distance = contacts[:,2].astype(int) - contacts[:,1].astype(int)
	distHist = np.histogram(distance, bins=(max(distance)))
	contact_counts = distHist[0]
	pos = np.arange(1, (max(distance))+1, 1)
	
	## Generate zero inflation model ##
	
	zeroPDF = []
	zero_prob = 0.0
	for k in range (0, DIST):
		zeros = 0
		if (k-50) < 0:
			zeroPDF.append(0.0)
		elif (k+50) > DIST:
				zeroPDF.append(zero_prob)
		else:
			for val in contact_counts[(k-50):(k+50)]:
				if val == 0:
					zeros += 1
			try:
				zero_prob = float(zeros)/100
			except:
				zero_prob = 0
			zeroPDF.append(zero_prob)

	zeroPDF = np.asarray(zeroPDF)
	zeroPos = np.arange(1, DIST+1, 1)
	lowess_sm = sm.nonparametric.lowess
	zeroSM = lowess_sm(zeroPDF,zeroPos,frac=0.01,it=3,delta=0,return_sorted = False)
	
	zeroModel = []
	for zeroProbability in zeroSM:
		zeroModel.append(0.5*(1-zeroProbability))
	zeroModel = np.asarray(zeroModel)
	
	## Generate short range model and append to bgModel ##
	#contact_counts = [x * 5 for x in contact_counts]#map(lambda x : x *.5, contact_counts)
	
	short_pos = pos[0:1000]
	short_counts = contact_counts[0:1000]
	pseudo_counts = [x + 0 for x in short_counts]#map(lambda x : x + 0, short_counts)
	short_pos = np.asarray(short_pos)
	pseudo_counts = np.asarray(pseudo_counts)
	lowess_sm = sm.nonparametric.lowess
	counts_sm = lowess_sm(pseudo_counts,short_pos,frac=0.05,it=3,delta=0,return_sorted = False)
	for i in counts_sm:
		bgModel.append(i)

	#Generate piecewise long range models ##
	for i in range(models):
		start,stop = ((i+1) * 50000)-50000, (i+1) * 50000
		if stop < DIST - 300:
			short_pos = pos[start:stop+300]
			short_counts = contact_counts[start:stop+300]
			zero_counts = zeroModel[start:stop+300]
			pseudo_counts = np.add(short_counts, zero_counts)
			short_pos = np.asarray(short_pos)
			pseudo_counts = np.asarray(pseudo_counts)
		else:
			short_pos = pos[start:stop]
			short_counts = contact_counts[start:stop]
			zero_counts = zeroModel[start:stop]
			pseudo_counts = np.add(short_counts, zero_counts)
			short_pos = np.asarray(short_pos)
			pseudo_counts = np.asarray(pseudo_counts)
		lowess_sm = sm.nonparametric.lowess
		counts_sm = lowess_sm(pseudo_counts,short_pos,frac=.01,it=3,delta=0,return_sorted = False)

	## Add long range model to bgModel ## 
		if i < 1:
			tail = bgModel[-300:]
			head = counts_sm[700:1001]
			merge = []
			for (x,y) in zip(tail,head):
				merge.append((x+y)/2.0)
			tmpModel,counts = bgModel[:-300],counts_sm[1001:]
			bgModel = []
			for val in tmpModel:   
				bgModel.append(val)
			for val in merge:
				bgModel.append(val)
			for val in counts:
				bgModel.append(val)
		else:
			tail = bgModel[-300:]
			head = counts_sm[:300]
			merge = []
			for (i,j) in zip(tail,head):
				merge.append((i+j)/2)
			tmpModel,counts = bgModel[:-300],counts_sm[300:]
			bgModel = []
			for val in tmpModel:   
				bgModel.append(val)
			for val in merge:
				bgModel.append(val)
			for val in counts:
				bgModel.append(val)
			
	bgPDF = []
	reads = sum(val < DIST for val in distance)

	print (reads)
	for base in bgModel:
		bgPDF.append(float(base/reads))
			
#		bgPDF = map(lambda x : x * 2, bgPDF)
#bgPDF = [x * 2 for x in bgPDF]

baitStart, baitStop = int(position-CAP), int(position+CAP)

## Get interactions around bait
for contact in contacts:
	if baitStart <= int(contact[1]) <= baitStop:
#					print ('plus hit!')
			plus.append(int(contact[2]))
	if baitStart <= int(contact[2]) <= baitStop:
			minus.append(int(contact[1]))
## test plus strand contacts with promoters
for i in preyFile_DF.index:
#		print ('prey0',  + str(prey[0]) + ' minDist: ' + str(minDist))
	if preyFile_DF['pos'][i] - position > minDist:
		preyChrom, preyPosition, preyDist = preyFile_DF['chrom'][i].split('chr')[1],preyFile_DF['pos'][i], preyFile_DF['pos'][i] - position
		preyStart, preyStop = (preyPosition - (CAP)), (preyPosition + (CAP))
		if (preyDist-CAP) >= 0:
			expStart, expStop = int(preyDist-CAP), int(preyDist + CAP)
			exp_prob = sum(bgPDF[expStart:(expStop+1)])
			expected = (len(plus)*exp_prob)
#				print ('Expected: ' + str(expected)) 
		else:
			expStart, expStop, mirror = 0, (preyDist+CAP), (-1*(preyDist-CAP))
			exp_prob = sum(bgPDF[expStart:expStop+1])
			exp_prob= exp_prob# + (sum(bgPDF[0:mirror+1]))
			expected = (len(plus)*exp_prob)
#				print ('Expected: ' + str(expected)) 
		
		observed = 0
		for connection in plus:
			if preyStart <= connection <= preyStop:
				observed += 1
#			print ('Plus: ' + str(len(plus)))
#			print str(position), str(len(plus)), str(observed)
			
		v = robjects.FloatVector([observed, expected, (len(plus)-observed), (len(plus)-expected)])
		robjects.r.assign("v", v)
		test = robjects.r("matrix(v, nrow = 2)")
		robjects.r.assign("matrix", test)
		try:
			robjects.r("exact <- fisher.test(matrix, alternative = 'greater')")
			result = robjects.r("exact")
			p_val = result[0][0]
#				print (result)
		except:
			p_val = 999.99
		contactProbabilities.append([chrom, position, preyChrom, preyPosition, preyDist, p_val, observed, expected])
			
## test minus strand contacts with promoters
for i in preyFile_DF.index:
	if preyFile_DF['pos'][i] - position < (-1*minDist):
		preyChrom, preyPosition, preyDist = preyFile_DF['chrom'][i].split('chr')[1],preyFile_DF['pos'][i], preyFile_DF['pos'][i] - position
		absDist = (preyDist * -1)
		preyStart, preyStop = (preyPosition - (CAP)), (preyPosition + (CAP))
		if (absDist-CAP) >= 0:
			expStart, expStop = (absDist-CAP), (absDist + CAP) 
			exp_prob = sum(bgPDF[expStart:expStop+1])
			expected = (len(minus)*exp_prob)
		else:
			expStart, expStop, mirror = 0, (absDist+CAP), (-1*(absDist-CAP))
			exp_prob = sum(bgPDF[expStart:expStop+1])
			exp_prob= exp_prob# + (sum(bgPDF[0:mirror+1]))
			expected = (len(minus)*exp_prob)
			
		observed = 0
		for connection in minus:
			if preyStart <= connection <= preyStop:
				observed += 1
#			print str(position), str(len(plus)), str(observed)
			
		v = robjects.FloatVector([observed, expected, (len(minus)-observed), (len(minus)-expected)])
		robjects.r.assign("v", v)
		test = robjects.r("matrix(v, nrow = 2)")
		robjects.r.assign("matrix", test)
		try:
			robjects.r("exact <- fisher.test(matrix, alternative = 'greater')")
			result = robjects.r("exact")
			p_val = result[0][0]
		except:
			p_val = 999.99
		contactProbabilities.append([chrom, position, preyChrom, preyPosition, preyDist, p_val, observed, expected])
		

out = open(sys.argv[4], "w")

print ('writing')
p_count = 0

for prob in contactProbabilities:
	# output file format - chr <\t> bait position <\t> prey position <\t> p_val <\t> observed <\t> expected <\t> FDR
	outStr = 'chr' + str(prob[0]) + "\t" + str(prob[1]) +"\t"+ str(prob[3]) +"\t"+ str(prob[4]) +"\t"+ str(prob[5]) +"\t"+ str(prob[6]) +"\t"+ str(prob[7]) +"\n"
	out.write(outStr)
	p_count += 1

