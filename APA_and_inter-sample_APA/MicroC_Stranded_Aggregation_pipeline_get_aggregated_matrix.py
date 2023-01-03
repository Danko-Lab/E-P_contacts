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
5 - heatmap output (svg)

'''

winSize = int(sys.argv[2])
pixNum = int(sys.argv[3])

matrix_list = os.listdir(sys.argv[1])

# starting with an all zero dataframe
aggregate = pd.read_csv(str(sys.argv[1]) + matrix_list[0], index_col = 0).astype(int)

# aggregate matrices
c = 0
for matrix in matrix_list[1:]:
	df = pd.read_csv(str(sys.argv[1]) + matrix, index_col = 0).astype(int)
	aggregate = sum([aggregate, df])
	c += 1
	print('adding matrices ' + str(c) + ' out of ' + str(len(matrix_list)))

#save aggregated matrix to a file
aggregate.to_csv(sys.argv[4])

#plot the aggregated matrix and save the APA heatmap
cmap = sns.color_palette("YlOrRd")
ax = sns.heatmap(aggregate, cmap = cmap, square = True)
plt.xticks([0,pixNum, pixNum*2], ['-' + str(int(winSize / 1000)) + 'kb', '0', str(int(winSize / 1000)) + 'kb'])
plt.yticks([0,pixNum, pixNum*2], [str(int(winSize / 1000)) + 'kb', '0', '-' + str(int(winSize / 1000)) + 'kb'])

plt.savefig(sys.argv[5])
