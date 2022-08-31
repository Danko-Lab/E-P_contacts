This directory includes the files required to calculate the ratio between contacts associated with each pair of enhancer
and promoter (or any two sets of anchors) in the genome, within a defined genomic distance range from each other, and their 
associated background contacts. The basis for this comperison is that the backgroung contacts are the contacts of the enhancer
with regions flanking the promoter and vice versa. 

![image](https://user-images.githubusercontent.com/47452349/187677535-e9bde4ff-0365-4bfa-8235-d375a74100b0.png)

Using this normalization strategy, we were able to ask whether enhancer-promoter or contacts between CTCF binding sites are 
lost between different states or conditions. We also include a python file where we visualized a comperison of these values 
for enhancer-promoter pairs between NELFB-FKBP12(F36V) mESCs treated or untreated with dTAG-13 for 30 minutes, using a scatterplot.

Usage:
    bash MicroC_EP_and_BG_contacts.bsh contacts_file baits_file preys_file minDis maxDis winSize OutputDir minBGdis maxBGdis
    
Required files and arguments:

contacts_file - contact (pairs) file output from distiler-nf filtered using:

    zcat perfix.pairs.gz | awk 'BEGIN {OFS = "\t"} ; {if ($1 == "." && $2 == $4 && $9 >= 30 && $10 >= 30) {print $2, $3, $4, $5, $6, $7, $8, $9, $10}}' > perfix.nodups_30_intra.pairs

baits_file - bed3 file with all the regions you want to use as baits (the script will center the heatmap over the center of the regions)

preys_file - bed3 file with all the regions you want to use as baits (the script will center the heatmap over the center of the regions)

minDis - the minimal distance between region centers to screen for. minDis needs to be >= winSize + max(frarment sizes in the library) (int)

maxDis - the maximal distance between region centers to screen for (int)

winSize - window half size of region to screen for contacts per region with other region in the list (int)

OutputDir - path to directory in which work will be done and output will appear

minBGdis - The minimal distance from the center of the second anchor to calculate background

maxBGdis - The maximal distance from the center of the second anchor to calculate background. 

NOTE: We made maxBGdis equal to maxDis in our study. Hence, the anchor will always be included in the range of its own background. However, we exclude the anchor itself from the background contacts (no intra anchor contacts) by limitations defined at the preys\baits_IN_contacts calculations.

Other requirements:

* The python file single_pair_contacts_and_background_calculation.py should be placed in the directory from which the bash script is being executed.

* Python>3

* Python packages:
  * sys
  * pandas
  * scipy
  * os
  * numpy
  * seaborn
  * matplotlib
