# Code for Barshad et al. 2023
This repository includes the main pieces of code used to support the conclusion in the paper. We also provide input files used in the paper when possible and explain how to obtain the ones that were too big to be provided here. The three main types of analyses are: Contact_normalization_by_local_decay, APA_and_inter-sample_APA and EP_contacts_compared_to_local_background. In each of the above repositories you can find a detailed explanation about the different scripts, what they do and requred\optional parameters and input files. Here, we provide examples for code execution with data analysed in this paper.   

Before we start the demo, please colne this repository and `cd` into it.

Dependencies:

- Python > 3
- Python packages:
  * sys
  * pandas
  * numpy
  * statsmodels
  * gzip
  * matplotlib
  * os
  * seaborn
  * scipy
- R packages:
  * rpy2

### Pairs files
In all of the provided code, the source for Micro-C contact information is from prosseced pairs files generated via the distiller-nf algorithm. Most of these processed files can be found in the GEO dataset for this project: GSE206131. Other are available at ftp://cbsuftp.tc.cornell.edu/danko/hub/MicroC_pairs_files/. A simple way of obtaining these files is to run distiller-nf with: `parsing_options: '--add-columns mapq' drop_readid: True`. You can find examples of raw Micro-C data processing in the "Micro-C_basic_processing" directory that containes YML configuration files used for distiller-nf. After obtaining pairs file for each replicate in your data, run:

    zcat perfix.rep1.pairs.gz perfix.rep2.pairs.gz perfix.rep3.pairs.gz ... | awk 'BEGIN {OFS = "\t"} ; {if ($1 == "." && $2 == $4 && $9 >= 30 && $10 >= 30) {print $2, $3, $4, $5, $6, $7, $8, $9, $10}}' > perfix.nodups_30_intra.pairs```

For the purpes of these demos, please make sure you have `cd` into the E-P_contacts repository and run the follwing to download the relevant processed pairs files:

    wget "ftp://cbsuftp.tc.cornell.edu/danko/hub/MicroC_pairs_files/*"
    wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206131/suppl/GSE206131_K562_cis_mapq30_pairs.txt.gz

### Contact_normalization_by_local_decay
We used this code to compare lacal background signal normalized contact frequency between functional and nonfunctional pairs of enhancers and promoters. The normalization to the local pattern of the distance decay provides a way of seperating the effect of contact frequency from that of genomic distance - two features that often highly correlate. To run this code, you will need the pairs file for our Micro-C data from K562 cells (GSE206131_K562_cis_mapq30_pairs.txt.gz) that can be found at the GEO repository. The rest of the input files can be found at the input files directory in this GitHub repository. To obtain observed and expected contacts between 4kb windows around enhancers and promoters that were tested by CRISPRi and are up to 1Mb from each other in K562 cells, run:

    bash ./Contact_normalization_by_local_decay/ContactCaller_microC.bsh ./Input_files/Gasperini_dREG_based_TRE_baits_hg38.txt ./Input_files/Gasperini_dREG_based_promoter_preys_hg38.txt GSE206131_K562_cis_mapq30_pairs.txt.gz ./Contact_normalization_by_local_decay/ 30 1000000 2000

Note: outputPath refers to the directory where work will be done. The acompaning python file, ContactCaller_microC.py, should be at the same directory as ContactCaller_microC.bsh. The above command will use 30 CPU cores.

After ContactCaller_microC.bsh finish running, run the following to concatenate the data for all enhancer-promoter pairs tested:

    cat ./Contact_normalization_by_local_decay/chr* > ./Contact_normalization_by_local_decay/ContactCaller_microC_output.txt

After getting the observed and expected contacts for each enhancer-promoter pair, we defined "functional", "nonfunctional" and "other" pairs, based on CRISPRi data. Ths can be done by:

    python ./Contact_normalization_by_local_decay/Pair_type_assignment.py ./Contact_normalization_by_local_decay/ContactCaller_microC_output.txt ./Input_files/Gasperini_dREG_based_functional.csv ./Input_files/Gasperini_dREG_based_nonfunctional.csv ./Contact_normalization_by_local_decay/ContactCaller_microC_output_W_functional_nonfunctional_and_other_pair_assignments.txt 

To then visualize the distribution of lacal decay-normalized contacts by pair type, while limiting for a minimum of 1 contact per pair and a minimum distance of 15kb, run:

    python ./Contact_normalization_by_local_decay/Plotting_obs_over_exp_distribution_by_pair_type.py ./Contact_normalization_by_local_decay/ContactCaller_microC_output_W_functional_nonfunctional_and_other_pair_assignments.txt 1 15000 ./Contact_normalization_by_local_decay/Violinplot_for_normalized_contacts_by_pair_type.svg

### APA_and_inter-sample_APA
This is an alternative approach to conduct an aggregated peak analysis (APA). It can be used both to plot the aggregated raw contacts between anchors in the data and to compare between the aggregated contacts in different samples (like control and treatment), while normalizing for changes in anchor-associated contacts (1D signal). Please refer to the README at the APA_and_inter-sample_APA or the methods section in the paper for a detailed explanation of why this is important in the case of transcriptional inhibitors and for a mathematical representation of the calculation included.

To run this code, you will need the pairs file for FLV, TRP or DMSO (control) treated mESCs from Hsieh et al., 2020 Mol. cell paper. We deposited a processed version of these pairs files in ftp://cbsuftp.tc.cornell.edu/danko/hub/MicroC_pairs_files/. You will also need the baits and prey files which can be found at the input files directory in this GitHub repository. To obtain the raw aggregated contacts 20kbX20kb matrix at 200bp resolution between enhancers and promoters within 25-150kb of genomic distance, as well as the vectors for the 1D signal around these enhancers and promoters at the DMSO control, use:

    bash ./APA_and_inter-sample_APA/MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh mESCs_DMSO_30_intra.mm10.nodups.pairs.gz ./Input_files/dREG_based_promoters_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed ./Input_files/dREG_based_TREs_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed 25000 150000 10000 50 ./APA_and_inter-sample_APA/DMSO

And for flavopiridol (FLV) treated cells: 

    bash ./APA_and_inter-sample_APA/MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh mESCs_FLV_30_intra.mm10.nodups.pairs.gz ./Input_files/dREG_based_promoters_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed ./Input_files/dREG_based_TREs_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed 25000 150000 10000 50 ./APA_and_inter-sample_APA/FLV

And for Triptolide (TRP) treated cells: 

    bash ./APA_and_inter-sample_APA/MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh mESCs_TRP_30_intra.mm10.nodups.pairs.gz ./Input_files/dREG_based_promoters_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed ./Input_files/dREG_based_TREs_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed 25000 150000 10000 50 ./APA_and_inter-sample_APA/TRP

After running MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh for all treatment and control samples, you can visualize the change APAs for 10530 (promoters) baits and 27900 (enhancers) preys as follows:

To get the flavopiridol (FLV) changes over DMSO control:

    python ./APA_and_inter-sample_APA/Change_calculation_and_visualization.py ./APA_and_inter-sample_APA/DMSO/AggMat.csv ./APA_and_inter-sample_APA/FLV/AggMat.csv ./APA_and_inter-sample_APA/DMSO/baits_genome_wide_contacts.csv ./APA_and_inter-sample_APA/DMSO/preys_genome_wide_contacts.csv ./APA_and_inter-sample_APA/FLV/baits_genome_wide_contacts.csv ./APA_and_inter-sample_APA/FLV/preys_genome_wide_contacts.csv 10000 50 10530 27900 ./APA_and_inter-sample_APA/FLV_over_DMSO_1D_normalized_change_APA.svg 

To get the triptolide (TRP) changes over DMSO control:

    python ./APA_and_inter-sample_APA/Change_calculation_and_visualization.py ./APA_and_inter-sample_APA/DMSO/AggMat.csv ./APA_and_inter-sample_APA/TRP/AggMat.csv ./APA_and_inter-sample_APA/DMSO/baits_genome_wide_contacts.csv ./APA_and_inter-sample_APA/DMSO/preys_genome_wide_contacts.csv ./APA_and_inter-sample_APA/TRP/baits_genome_wide_contacts.csv ./APA_and_inter-sample_APA/TRP/preys_genome_wide_contacts.csv 10000 50 10530 27900 ./APA_and_inter-sample_APA/TRP_over_DMSO_1D_normalized_change_APA.svg 

### EP_contacts_compared_to_local_background
This is a complementary method to the APA, in which rather than examining the aggregated contacts, we are looking at the distribution of contacts across all pairs of enhancers and promoter (or other sets of enchors) within a defined range of genomic distances, while normalizing for the contacts obtained between the enhancer and the promoter nearby regions and vise-versa. The benifit of this approach is that it is, on one hand captures the entire distribution of changes between control and treatment conditions, while on the other hand, less affected by outliers. It also allows us to perform ststistical tests to ask if the global trend obtained under different treatment is statistically significant.

To run this code, you will need the pairs file for FLV, TRP or DMSO (control) treated mESCs from Hsieh et al., 2020 Mol. cell paper. We deposited a processed version of these pairs files in ftp://cbsuftp.tc.cornell.edu/danko/hub/MicroC_pairs_files/. You will also need the baits and prey files which can be found at the input files directory in this GitHub repository. To obtain the contacts between 5kb windows around enhancers and promoters within 25-150kb of genomic distance and background regions being 10-150kb away from enhancers and promoters TSS, at the DMSO control, run:

    bash ./EP_contacts_compared_to_local_background/MicroC_EP_and_BG_contacts.bsh mESCs_DMSO_30_intra.mm10.nodups.pairs.gz ./Input_files/dREG_based_promoters_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed ./Input_files/dREG_based_TREs_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed 25000 150000 10000 150000 ./EP_contacts_compared_to_local_background/DMSO/
    
And for flavopiridol (FLV) treated cells: 

    bash ./EP_contacts_compared_to_local_background/MicroC_EP_and_BG_contacts.bsh mESCs_FLV_30_intra.mm10.nodups.pairs.gz ./Input_files/dREG_based_promoters_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed ./Input_files/dREG_based_TREs_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed 25000 150000 10000 150000 ./EP_contacts_compared_to_local_background/FLV/

And for triptolide (TRP) treated cells: 

    bash ./EP_contacts_compared_to_local_background/MicroC_EP_and_BG_contacts.bsh mESCs_TRP_30_intra.mm10.nodups.pairs.gz ./Input_files/dREG_based_promoters_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed ./Input_files/dREG_based_TREs_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed 25000 150000 10000 150000 ./EP_contacts_compared_to_local_background/TRP/

Note: outputPath refers to the directory where work will be done. The acompaning python file, single_pair_contacts_and_background_calculation.py should be at the same directory as MicroC_EP_and_BG_contacts.bsh.

After running MicroC_EP_and_BG_contacts.bsh for all three treatmen and control conditions, you can run the following to obtain scatterplots comparing the EP ovver background ratios across all EP pairs, between the different treatment conditions:

    python ./EP_contacts_compared_to_local_background/Compering_EP_contacts_between_treatments.py ./EP_contacts_compared_to_local_background/DMSO/EP_and_BG_contacts.txt ./EP_contacts_compared_to_local_background/FLV/EP_and_BG_contacts.txt ./EP_contacts_compared_to_local_background/TRP/EP_and_BG_contacts.txt 53226768 362862200 410040533 8 ./EP_contacts_compared_to_local_background/FLV_vs_DMSO.svg ./EP_contacts_compared_to_local_background/TRP_vs_DMSO.svg ./EP_contacts_compared_to_local_background/FLV_vs_TRP.svg

Note: The integers above represent (from left to righr) the sequencing depth of DMSO, FLV and TRP treated mESCs Micro-C for cis interactions with mapq>30 and the minimum threshold for contacts per billion (CPB) for EP contact, to include in the analysis.
