# This script will sum all individual bait's APAs for the final, genome-wide APA for the sample

import matplotlib
matplotlib.use('SVG')
import sys
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt

'''
1 - path to matrices directory
2 - window half size of region to screen for contacts per region with other region in the list (integer)
3 - half the number of pixels to devide the window to (integer). The resulting pixels will represent a window-size/number-of-pixels bps
4 - aggregated matrix output (csv)

'''

winSize = int(sys.argv[2])
pixNum = int(sys.argv[3])


matrix_list = os.listdir(sys.argv[1])

# starting with an all zero dataframe
aggregate = pd.read_csv(str(sys.argv[1]) + matrix_list[0], index_col = 0).astype(int)

c = 0
for matrix in matrix_list[1:]:
	df = pd.read_csv(str(sys.argv[1]) + matrix, index_col = 0).astype(int)
	aggregate = sum([aggregate, df])
	c += 1
	print('adding matrices ' + str(c) + ' out of ' + str(len(matrix_list)))

aggregate.to_csv(sys.argv[4])

