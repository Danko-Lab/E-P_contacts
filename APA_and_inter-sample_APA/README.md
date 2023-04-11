*Overview*

We expected significant changes in chromatin after manipulating Pol II transcription. As such, not only are enhancer-promoter contacts expected to change, but the background contacts with at least one end originating at enhancer- and promoter- regions may be affected between conditions. As APAs are often used to characterize contacts, we devised an APA that normalizes enhancer-promoter contacts for changes in the 1D signal mapping to either anchor region. The primary challenge with devising a background-corrected APA is to handle the sparsity of Micro-C data (i.e., most small genomic bins have an observed contact value of 0). To address this challenge, our strategy computes a single observed and background matrix separately that represents the set of all enhancer and promoter regions included in the analysis, and then performs a pixel-by-pixel division of the aggregate observed and background matrices  (Fig. 2A, Fig. S5). 

*Computing the observed matrix*

We first compute an aggregate observed matrix that represents all enhancer-promoter pairs.  To do this calculation, we take the sum over the set of all enhancer-promoter pair matrices, leaving a single aggregate observed matrix (represented in Fig. S5A). For a single enhancer-promoter pair, , we calculated an observed contact matrix, $\boldsymbol{O}^{k}$: 

$$\boldsymbol{O}^{k} = [C^{k}_{i,j}]$$ 

Where $C^{k}_{i,j}$ is the number of contacts mapped to the $i$th window relative to the enhancer TSS and the $j$th window relative to the promoter TSS for enhancer-promoter pair $k$. 

To make a single aggregate matrix, we take the element-wise sum of the matrices for all enhancer-promoter pairs. The result is a single observed ($\boldsymbol{obs}$) APA matrix that represents the aggregate signal across all enhancer-promoter pairs (Fig. S5A). Formally, the computation is was completed as follows:

$$\boldsymbol{obs} = \sum_{k=1}^{n_{ep}}\boldsymbol{O}^{k}$$

Where $n_{ep}$, the number of enhancer-promoter pairs, is restricted by the allowed genomic distance range we defined. In figures calculating APAs within a single sample or condition, we presented the $\boldsymbol{obs}$ matrix as heatmaps (Figs. 1D and 3B). 

*Computing the background matrix*

As noted above, different conditions can have very different densities of Micro-C paired-end tags mapping to each anchor, which are in some cases caused by effects the perturbation studied in the experiment has on nucleosome occupancy and/ or other aspects of 1D chromatin structure. Interpreting whether enhancer-promoter interactions change between conditions requires asking whether the observed change is more extreme than can be explained by changes in the density of 1D background signal mapping at both anchors (henceforth called background). Therefore when representing the changes in contacts between different conditions, we normalized the observed changes to background changes in 1D signal near enhancers and promoters.

We computed the background matrix in two steps: (1) We compute the 1D signal near each enhancer/ promoter anchor, and (2) We turn the 1D signal into a matrix by computing the outer sum of the signal at enhancer/ promoter anchors. The motivation for this strategy is that we assume the probability of observing a signal in window $i,J$ of the matrix based on the 1D signal is proportional to the probability of observing a read in either window $i$ in the enhancer or window $j$ in the promoter. These assumptions motivate the use of the sum of signals in each anchor in each window to build the matrix (often called the outer sum).

Next, we compute a background matrix that represents the aggregate signal near the enhancer-promoter anchors in the same dataset. We compute the background matrix in two steps: (1) We compute the average 1D signal near each enhancer and each promoter anchor, and (2) We turn the average 1D signal into a matrix by computing the outer sum of the signal at enhancer and promoter anchors.

The motivation for this strategy is that we assume the probability of observing a signal in window $i,j$ of the observed matrix is proportional to the probability of observing a read in either window $i$ in the enhancer or window $j$ in the promoter. Further we assume the probabilities observing reads in window $i$ and $j$ are statistically independent. These assumptions motivate the use of the sum of signals in each anchor in each window to build the matrix, often called the outer sum, because the probability of observing a read from either the enhancer or promoter is the sum of the two probabilities. We also considered alternative formulations that convert 1D signal vectors into a matrix using the outer product, instead of the outer sum. The problem with this formulation in the setup used here is that the outer product includes terms for potential enhancer-promoter pairs that do not meet the criteria used in our analysis, and therefore were not incorporated in the observed (e.g., including cases where the enhancer-promoter pair reside on different chromosomes) (see Supplementary Note 1).
 
The computation of the background matrix is performed as follows: 

First (step 1), we defined a vector of counts with the same length/ width as the APA matrix, $\boldsymbol{O}$, that represents the sum of all Micro-C paired-end tags in which at least one end falls into that window relative to the anchor (usually the enhancer or promoter TSS), and the other end falls between the minimum and maximum distance allowed between enhancers and promoters in the APA. Formally, we first compute vectors that represent the aggregate 1D signal at positions $i$ or $j$ for enhancers or promoters, respectively. We take the mean signal over the set of all $n_{e}$ enhancers or $n_{p}$ promoters in the dataset, as shown:

$$E_{i} = \sum_{x=1}^{n_{e}}E^{x}_{i}/n_{e}$$

$$P_{j} = \sum_{y=1}^{n_{p}}P^{y}_{j}/n_{p}$$

These vectors are shown in the bottom right panel of Fig. S5B. Note that we index using $x$ and $y $(instead of $k$) to emphasize that these reflect individual enhancers and promoters, rather than enhancer-promoter pairs.

Second (step 2), we used vectors $E$ and $P$  to generate the background matrix $\boldsymbol{B}$. We compute the background matrix using the outer sum. Hence, the calculation of cell $i,j$ in the background matrix is computed as follows: 

$$\boldsymbol{B}_{i,j} = E_{i} + P_{j}$$

Or, in vector notation, matrix $\boldsymbol{B}$ is defined as the outer sum of vectors $E$ and $P$:

$$\boldsymbol{B} = E \bigoplus P$$

Note that $n_{ep}$ is not related to either $n_{e}$ or $n_{p}$ for two reasons: first, not all enhancer-promoter pairs are allowed by our distance requirements, and second each enhancer (or promoter) can be paired with multiple promoters (or enhancers). 

*Computing background corrected APAs*

We calculated the background corrected APA matrix $\boldsymbol{F}$ by dividing observed and background matrices for each condition, as follows: 

$$\boldsymbol{F}_{i,j} = \frac{\boldsymbol{obs}_{i,j} (T) / \boldsymbol{B}_{i,j} (T)}{\boldsymbol{obs}_{i,j} (C) / \boldsymbol{B}_{i,j}(C)}$$

Where $T$ stands for treatment condition and $C$ for control condition. 

*Heuristics for choosing which enhancer-promoter pairs are included in each analysis*

The primary concern when choosing window sizes in the APA is to avoid overlapping windows between enhancers and promoters, which would result in crossing the diagonal of the Micro-C matrix. To avoid overlapping windows around enhancers and promoters, we excluded enhancer-promoter pairs for which the separating genomic distance was smaller than the total 1D size of the APA plus the maximal fragment size in the library, which is known for Micro-C libraries due to the agarose gel purification step. For example, for a 20kb x 20kb APA, the minimum enhancer-promoter distance should be larger than 20.3 kb. For APAs calculated at windows of 20kb around the anchors, we considered all possible anchor pairs within a genomic distance of 25-150kb. For the high resolution APA with 2kb window around enhancer and promoter TSSs (Fig. 1D), we considered all possible enhancer-promoter pairs within a genomic distance of 5-100kb.


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

It should be ran as following:

    python Change_calculation_and_visualization.py control_raw_APA treatment_raw_APA control_baits_1D_signal control_preys_1D_signal treatment_baits_1D_signal treatment_preys_1D_signal winSize pixNum number_of_baits number_of_preys path_to_APA_image_file

Required files and arguments:

control_raw_APA - the control APA matrix (a csv file output of MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh)

treatment_raw_APA - the control APA matrix (a csv file output of MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh)

control_baits_1D_signal - the control 1D signal vector around baits (a csv file output of MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh)

control_preys_1D_signal - the control 1D signal vector around preys (a csv file output of MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh)

treatment_baits_1D_signal - the treatment 1D signal vector around baits (a csv file output of MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh)

treatment_preys_1D_signal - the treatment 1D signal vector around preys (a csv file output of MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh)

winSize - window half size of region to screen for contacts per region with other region in the list (integer). This will be half the APA 1D size.

pixNum - half the number of pixels to devide the window to (integer). The resulting pixels will represent winSize/pixNum bps

number_of_baits - an integer representing the number of baits in the data

number_of_preys - an integer representing the number of preys in the data

path_to_APA_image_file

Other requirements:

* The python files MicroC_Stranded_Aggregation_pipeline_get_bait_matrix.py, MicroC_Stranded_Aggregation_pipeline_get_aggregated_matrix.py and get_genome_wide_normalization_scores_by_search_window should be placed in the directory from which the bash script MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh is being executed.
