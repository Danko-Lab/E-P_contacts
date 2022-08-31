This directory holds the bash and associated python files to calculate the intra-sample aggregated peak analysis (APA),
calculate the 1D signal associated with each set of anchors (baits and preys) in each sample and calculate the ratio
of ratios representing observed over expected change in aggregated contacts between the anchors.

Usage:

The first step is to calculate the APA for each sample, using the contacts file output of distiller-nf. To focus on
nucleosome contacts, thi code centers the reads in the contacts file by shifting them 75bp inward as such:

 ![image](https://user-images.githubusercontent.com/47452349/187545681-095b5e4c-fc92-4491-9ec7-06a3790b757e.png)

Then, the APA is calculated using the following command:

    bash MicroC_Stranded_Aggregation_Pipeline.bsh contacts_file baits_file preys_file minDis maxDis winSize pixNum OutputDir

Required files and arguments:

contacts_file - contact (pairs) file output from distiler-nf filtered using:

    zcat perfix.pairs.gz | awk 'BEGIN {OFS = "\t"} ; {if ($1 == "." && $2 == $4 && $9 >= 30 && $10 >= 30) {print $2, $3, $4, $5, $6, $7, $8, $9, $10}}' > perfix.nodups_30_intra.pairs

baits_file - bed3 file with all the regions you want to use as baits (the script will center the heatmap over the center of the regions)

preys_file - bed3 file with all the regions you want to use as baits (the script will center the heatmap over the center of the regions)

minDis - range start: the minimal distance between region centers to screen for. minDis needs to be >= winSize + max(frarment sizes in the library)

maxDis -  range end: the maximal distance between region centers to screen for

winSize - window half size of region to screen for contacts per region with other region in the list (integer). This will be the APA 1D size.

pixNu - half the number of pixels to devide the window to (integer). The resulting pixels will represent winSize/pixNum bps

OutputDir - path to directory in which work will be done and output will appear

The second step is also executed by a bash script and calculates the the aggregated 1D contact signal associated with each anchor set.
This step will need to be done for each of the anchor sets seperatly, using the following command:

    bash get_genome_wide_normalization_scores_for_inter-sample_pileup_analysis.bsh anchor_type winSize pixNum OutputDir
    
Required files and arguments:

anchor_type - "baits" of "prays" (str)

winSize - same as winSize in step 1

pixNum - same as pixNum in step 1

OutputDir - same as OutputDir in step 1 (This is crucial for the code to run without errors!)

The third step in calculating an inter-samples contact change APA is to calculate the expected APA matrix based on the above calculated 1D signal,
devide the calculate a cell-by-cell (pixel-py-pixel) ratio between the observed and between the expected matrices of the two samples, and generate
the observed/expected change APA matrix by cell-by-cell dividing the observed and expected change matrices. An example of how it was done for the 
comperison of the flavopiradol- or triptolide-treated mESCs to the DMSO control taken from Hsieh et al. 2020 (GSE130275) is provided in the python
file named Change_calculation_and_visualization.py. this file includes an example for the APA obs/exp change matrices calculation based on the outputs
of steps 1 and 2 above, as well as visualization of the APA matices as heatmaps, smoothing of these matrices and lotting a boxplot with the calculated
values associated with the dot (center), the stripes stemming from the dot and the edges of the APA matrix. A clear description of every step is provided
in the python file.


Other requirements:

* The python files MicroC_Stranded_Aggregation_pipeline_get_bait_matrix.py and MicroC_Stranded_Aggregation_pipeline_get_aggregated_matrix.py should be placed in the directory from which the bash script MicroC_Stranded_Aggregation_pipeline.bsh is being executed.

* The python files get_genome_wide_normalization_scores_for_inter-sample_pileup_analysis_single_chromosome.py and get_genome_wide_normalization_scores_for_inter-sample_pileup_analysis_combine_chromosomes.py should be placed in the directory from which the bash script get_genome_wide_normalization_scores_for_inter-sample_pileup_analysis_single_chromosome.bsh is being executed.

* Python>3

* Python packages:
  * sys
  * pandas
  * matplotlib
  * os
  * seaborn
  * numpy
  * scipy
