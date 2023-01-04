This directory includes the files required to calculate the ratio between contacts associated with each pair of enhancer
and promoter (or any two sets of anchors) in the genome, within a defined genomic distance range from each other, and their 
associated background contacts. 

This is an alternative normalization scheme to our inter-sample APA, which compares the number of contacts between enhancer-promoter pairs to the local background near each enhancer and promoter anchor. The primary goal of this alternative normalization scheme was to assess changes in contact frequency specific to enhancer-promoter pairs, after accounting for changes in contact frequency between the enhancer (or promoter) and flanking regions to the second anchor. We calculated the number of contacts between each pair of anchors (enhancer-promoter or CTCF binding sites) using a 5kb window around each anchor. As a background, we counted  the number of contacts between each anchor (in a 5kb window) and regions 10-150 kb from the second anchor (Fig. 2C). As such, the ratio between the anchor-to-anchor (i.e., either enhancer-promoter or CTCF-CTCF) contacts and background contacts was calculated for each pair of anchors using the following formula:

$$Background \text{ } normalized  \text{ } contacts = \frac{E \in P}{E \in Pbg + P \in Ebg}$$

Where $E \in P$ represents all contacts mapped to a 5kb window around both the enhancer and promoter TSSs, respectively. $E \in Pbg$ and $P \in Ebg$ represent contacts between 5k of the enhancer (or promoter) and the background window (10-150bp) of the promoter (or enhancer), respectively.

These background normalized contacts are computed separately for each anchor pair and are presented in scatterplots, box-and-whiskers plots or line plots over the NELFB degradation and dTAG washout time course, to calculate the statistical significance of changes between treatments and samples. To avoid the impact of noise, we analyzed only contacts that met a minimum baseline of anchor-to-anchor contacts (at least 8 contacts per billion contacts (CPB)) in one of the treatment conditions. Since TSS calling data (PRO-seq and coPRO-capped) was more abundant for K562 than Jurkat, when comparing K562 and Jurkat libraries we considered enhancer-promoter pairs with at least 8 CPB in both cell lines, to avoid ascertainment bias. The distribution of ratios between enhancer-promoter and background contacts in treated samples (Olaparib\TRP\FLV or dTAG treated cells) was compared to the median ratio in the respective control samples (Figs 2D-E, 4C-D, 5B and S6A). For comparison between cell lines, enhancer-promoter contacts at promoters with increased gene body transcription and\or Pol II pausing signal in one cell line, were compared to their median at the other cell line (Figs. 3D-E and S4A,C).

To avoid overlap, we had set the minimum genomic distance between anchors in each pair to 25kb. Additionally, as many anchors can be included in the background-associated regions flanking the second anchor, we excluded any contacts where both ends fall within the anchor’s defined window (intra-anchor contact) from the background contacts.   


![image](https://user-images.githubusercontent.com/47452349/187677535-e9bde4ff-0365-4bfa-8235-d375a74100b0.png)

Using this normalization strategy, we were able to ask whether enhancer-promoter or contacts between CTCF binding sites are 
lost between different states or conditions. We also include a python file where we visualized a comperison of these values 
for enhancer-promoter pairs between $\text{NELFB-FKBP12}^{F36V}$ mESCs treated or untreated with dTAG-13 for 30 minutes, using a scatterplot.

Usage:
    bash MicroC_EP_and_BG_contacts.bsh contacts_file baits_file preys_file minDis maxDis winSize OutputDir minBGdis maxBGdis
    
Required files and arguments:

contacts_file - contact (pairs) file output from distiler-nf without pair name and with mapq column, which is filtered to keep only cis contacts with both sides having a mapq >30 using:

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
