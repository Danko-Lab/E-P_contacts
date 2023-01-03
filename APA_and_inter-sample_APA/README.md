This directory holds the bash and associated python files to calculate the intra-sample aggregated peak analysis (APA),
calculate the 1D signal associated with each set of anchors (baits and preys) in each sample and calculate the ratio
of ratios representing observed over expected change in aggregated contacts between the anchors.

We expected significant changes in chromatin after manipulating Pol II transcription. As such, not only are enhancer-promoter contacts expected to change, but the background contacts with at least one end originating at enhancer- and promoter- regions may be affected between conditions. As APAs are often used to characterize contacts, we devised an APA that normalizes enhancer-promoter contacts to all contacts associated with enhancers and promoters.

Our strategy normalized contacts in the APA by changes in aggregated 1D signal mapped to the same windows used in the APA between conditions. We calculated the expected change in each pixel based on the ratio between the sum of the 1D signal in each sample or treatment condition (Fig. 2A). We computed the enhancer-promoter contact in each  in each pixel relative to the TSS of enhancers and promoters as:

$$obs_{i,j} = \sum_{k=1}^{n_{ep}} M_{(k)i',j'}$$

Where the signal at pixel $i$ (relative to the enhancer TSS) and $j$ (relative to the promoter TSS) is summed across the $i’$ and $j’$ pixels in each of the n matrices ( $M$ ) representing enhancer-promoter pairs within the defined genomic distance range from each other. Each matrix $(M_{(k)})$ represents the contact matrix between a single enhancer-promoter pair.

The background expected signal in each pixel, for each sample, was calculated as the 1D signal density relative to the enhancer and promoter TSS. This quantity was calculated using the following formula:

$$exp_{i,j} = \sum_{k=1}^{n_{e}}E_{(k)i'} / n_e + \sum_{k=1}^{n_{p}}P_{(k)j'} / n_p$$

Where $n_e$ is the total number of enhancers in the data and $n_p$ is the total number of promoters in the data. The 1D signal density at position $i$ relative to the promoter and $j$ relative to the enhancer ( $E(k)i'$ ) and each  promoter ( $P(k)j'$ ) was read and summed across all enhancers and promoters. The minimal distance of contacts considered at the 1D signal calculations was the same as in the enhancer-promoter contacts calculation, to avoid noise stemming from random ligation events for loci at very close linear proximity. Note that we have divided the sums over position $i$ relative to the promoter and $j$ relative to the enhancer by $n_e$ and $n_p$, respectively to avoid a bias due to a higher number of enhancers in the data. 


To compare the 1D normalized changes in contacts between treatment conditions, we calculated a ratio of ratios, with the numerator ratio representing the observed change and the denominator representing the expected change in contacts for that pixel. So, to account for the obs/exp changes in samples A relative to sample B, we used the following formula:

$$obs / exp = \frac{obsA / obsB}{expA / expB}$$

For all other, intra-sample APAs we used the aggregated raw counts, as described for the observed values calculation above. 

To avoid overlapping windows around enhancers and promoters, we excluded enhancer-promoter pairs for which the separating genomic distance was smaller than the total 1D size of the APA plus the maximal fragment size in the library, which is known for Micro-C libraries due to the agarose gel purification step. For example, for a 20kb x 20kb APA, the minimum enhancer-promoter distance should be larger than ~20.3 kp. For APAs calculated at windows of 20kb around the anchors, we considered all possible anchor pairs within a genomic distance of 25-150kb. For higher resolution APAs with 2kb window around enhancer and promoter TSSs (Fig. 1D and Fig. S3), we considered all possible enhancer-promoter pairs within a genomic distance of 5-100kb.



Usage:

The bash script MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh takes a pairs file as descibed below, baits and preys files, as well as numerical parameters for the minimum and maximum genomic distance between baits and preys, the 1D size of the APA matrix and the number of cells (or pixels) in each row and column and outpudts the raw contact APA matrix as well as vectors for the 1D signal around baits and preys.  

It should be ran as following:

    bash MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh contacts_file baits_file preys_file minDis maxDis winSize pixNum OutputDir

Required files and arguments:

contacts_file - contact (pairs) file output from distiler-nf without pair name and with mapq column, which is filtered to keep only cis contacts with both sides having a mapq >30 using:

    zcat perfix.pairs.gz | awk 'BEGIN {OFS = "\t"} ; {if ($1 == "." && $2 == $4 && $9 >= 30 && $10 >= 30) {print $2, $3, $4, $5, $6, $7, $8, $9, $10}}' > perfix.nodups_30_intra.pairs

baits_file - bed3 file with all the regions you want to use as baits (the script will center the heatmap over the center of the regions)

preys_file - bed3 file with all the regions you want to use as baits (the script will center the heatmap over the center of the regions)

minDis - range start: the minimal distance between region centers to screen for. minDis needs to be >= winSize + max(frarment sizes in the library)

maxDis -  range end: the maximal distance between region centers to screen for

winSize - window half size of region to screen for contacts per region with other region in the list (integer). This will be half the APA 1D size.

pixNu - half the number of pixels to devide the window to (integer). The resulting pixels will represent winSize/pixNum bps

OutputDir - path to directory in which work will be done and output will appear

Note: To focus on nucleosome contacts, the code centers the reads in the contacts file by shifting them 75bp inward as such:

 ![image](https://user-images.githubusercontent.com/47452349/187545681-095b5e4c-fc92-4491-9ec7-06a3790b757e.png)

Finally, the python script Change_calculation_and_visualization.py takes the raw APA matrices and 1D signal vectors around baits and preys from a treatment and control datasets, as well as numerical parameters for the 1D size of the APA matrix and the number of cells (or pixels) in each row and column and outpudts the matrix of 1D signal normalized change APA for contacts between the baits and preys.    

Other requirements:

* The python files MicroC_Stranded_Aggregation_pipeline_get_bait_matrix.py, MicroC_Stranded_Aggregation_pipeline_get_aggregated_matrix.py and get_genome_wide_normalization_scores_by_search_window should be placed in the directory from which the bash script MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh is being executed.

* Python>3

* Python packages:
  * sys
  * pandas
  * matplotlib
  * os
  * seaborn
  * numpy
  * scipy
