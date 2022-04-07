### This script is for combining the total contacts mapped to each of the defined windows
### aroung baits\preys that were calculated for single chromosomes


import sys
import pandas as pd

'''
1 - output (and input) directory
'''
#start with the chr1 table
aggregate = pd.read_csv(str(sys.argv[1]) + 'per_chromosome_normalization_tables/chr1.csv')
#aggregate all the rest of the chromosome (chromsome names can be changes for different species)
for chr in ['chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']:
	aggregate = sum(aggregate, pd.read_csv(str(sys.argv[1]) + 'per_chromosome_normalization_tables/' + chr + '.csv'))
#save the aggregated table of global contacts mapped to windows around baits\preys genome-wide
aggregate.to_csv(str(sys.argv[1]) + 'genome_wide_normalization_scores_table.csv')

