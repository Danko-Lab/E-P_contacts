This directory holds the bash and associated python files to calculate the intra-sample aggregated peak analysis (APA),
calculate the 1D signal associated with each set of anchors (baits and preys) in each sample and calculate the ratio
of ratios representing observed over expected change in aggregated contacts between the anchors.

We expected significant changes in chromatin after manipulating Pol II transcription. As such, not only are enhancer-promoter contacts expected to change, but the background contacts with at least one end originating at enhancer- and promoter- regions may be affected between conditions. As APAs are often used to characterize contacts, we devised an APA that normalizes enhancer-promoter contacts to all contacts associated with enhancers and promoters.

Our strategy normalized contacts in the APA by changes in aggregated 1D signal mapped to the same windows used in the APA between conditions. We calculated the expected change in each pixel based on the ratio between the sample-specific sum of the 1D signal in each sample or treatment condition (Fig. 2A). We computed the changes in 1D signal in each pixel in the APA as:

$$obs = \sum_{k=1}^n E_{i,k} \in P_{j,k}$$

Where $n$ is the total number of enhancer-promoter pairs, $E_{i,k} \in P_{j,k}$ are the contacts that fall within pixel $i$ (relative to the enhancer TSS) and $j$ (relative to the promoter TSS) in pair $k$. 

The expected contact values in each pixel, for each sample, were calculated as the expected read density based on 1D signal relative to the enhancer and promoter TSS. This quantity was calculated using the following formula:

$$exp = \left( \left( \sum_{k=1}^{n_e} E_{i,k} \right) / n_e \right) + \left( \left( \sum_{k=1}^{n_p} P_{j,k} / n_p \right) \right)$$

Where $n_e$ is the total number of enhancers in the data and $n_p$ is the total number of promoters in the data. Note that for the expected values calculation we have divided the sums over $E_i$ and $P_j$ by the total number of enhancers and promoters, respectively, to avoid a bias due to a higher number of enhancers in the data.

We then, for each pixel, calculated a ratio of ratios, with the numerator ratio representing the observed change and the denominator representing the expected change in contacts for that pixel. So, to account for the obs/exp changes in samples A relative to sample B, we used the following formula:

$$obs / exp = \frac{obsA / obsB}{expA / expB}$$

For all other, intra-sample APAs we used the aggregated raw counts, as described for the observed values calculation above. 

To avoid overlapping windows around enhancers and promoters, we excluded enhancer-promoter pairs for which the distance between them was smaller than the total 1D size of the APA plus the maximal fragment size in the library, which is known for Micro-C libraries due to the agarose gel purification step. For example, for a 20kb x 20kb APA, the minimum enhancer-promoter distance should be larger than ~20.3 kp. For APAs calculated at windows of 20kb around the anchors, we considered all possible anchor pairs within a genomic distance of 25-150kb. For higher resolution APAs with 2kb window around enhancer and promoter TSSs (Fig. 1D and Fig. S3), we considered all possible enhancer-promoter pairs within a genomic distance of 5-100kb.


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
