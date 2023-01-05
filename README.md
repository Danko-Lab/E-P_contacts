# Code for Barshad et al. 2023
This repository includes the main pieces of code used to support the conclusion in the paper. We also provide input files used in the paper when possible and explain how to obtain the ones that were too big to be provided here. The three main types of analyses are: Contact_normalization_by_local_decay, APA_and_inter-sample_APA and EP_contacts_compared_to_local_background. In each of the above repositories you can find a detailed explanation about the different scripts, what they do and requred\optional parameters and input files. Here, we provide examples for code execution with data analysed in this paper.   

### Pairs files
In all of the provided code, the source for Micro-C contact information is from prosseced pairs files generated via the distiller-nf algorithm. Most of these processed files can be found in the GEO dataset for this project: GSE206131. Other are available at ftp://cbsuftp.tc.cornell.edu/danko/hub/MicroC_pairs_files/. A simple way of obtaining these files is to run distiller-nf with: `parsing_options: '--add-columns mapq' drop_readid: True`. You can find examples of raw Micro-C data processing in the "Micro-C_basic_processing" directory that containes YML configuration files used for distiller-nf. After obtaining pairs file for each replicate in your data, run:

``zcat perfix.rep1.pairs.gz perfix.rep2.pairs.gz perfix.rep3.pairs.gz ... | awk 'BEGIN {OFS = "\t"} ; {if ($1 == "." && $2 == $4 && $9 >= 30 && $10 >= 30) {print $2, $3, $4, $5, $6, $7, $8, $9, $10}}' > perfix.nodups_30_intra.pairs```

### Contact_normalization_by_local_decay
We used this code to compare lacal background signal normalized contact frequency between functional and nonfunctional pairs of enhancers and promoters. The normalization to the local pattern of the distance decay provides a way of seperating the effect of contact frequency from that of genomic distance - two features that often highly correlate. To run this code, you will need the pairs file for our Micro-C data from K562 cells (GSE206131_K562_cis_mapq30_pairs.txt.gz) that can be found at the GEO repository. The rest of the input files can be found at the input files directory in this GitHub repository. To obtain observed and expected contacts between 4kb windows around enhancers and promoters that were tested by CRISPRi and are up to 1Mb from each other in K562 cells, run:

`bash ContactCaller_microC.bsh Gasperini_dREG_based_TRE_baits_hg38.txt Gasperini_dREG_based_promoter_preys_hg38.txt GSE206131_K562_cis_mapq30_pairs.txt.gz outputPath 30 1000000 2000`

Note: outputPath refers to the directory where work will be done. The acompaning python file, ContactCaller_microC.py, should be at the same directory as ContactCaller_microC.bsh. The above command will use 30 CPU cores.

After getting the observed and expected contacts for each enhancer-promoter pair, we defined "functional", "nonfunctional" and "other" pairs, based on CRISPRi data. Ths can be done by:

`python Pair_type_assignment.py ContactCaller_microC_output.txt Gasperini_dREG_based_functional.csv Gasperini_dREG_based_nonfunctional.csv ContactCaller_microC_output_W_functional_nonfunctional_and_other_pair_assignments.txt` 

To then visualize the distribution of lacal decay-normalized contacts by pair type, while limiting for a minimum of 1 contact per pair and a minimum distance of 15kb, run:

`python Plotting_obs_over_exp_distribution_by_pair_type.py ContactCaller_microC_output_W_functional_nonfunctional_and_other_pair_assignments.txt 1 15000 Violinplot_for_normalized_contacts_by_pair_type.svg`

###APA_and_inter-sample_APA
This is an alternative approach to conduct an aggregated peak analysis (APA). It can be used both to plot the aggregated raw contacts between anchors in the data and to compare between the aggregated contacts in different samples (like control and treatment), while normalizing for changes in anchor-associated contacts (1D signal). Please refer to the README at the APA_and_inter-sample_APA or the methods section in the paper for a detailed explanation of why this is important in the case of transcriptional inhibitors and for a mathematical representation of the calculation included.

To run this code, you will need the pairs file for FLV, TRP or DMSO (control) treated mESCs from Hsieh et al., 2020 Mol. cell paper. We deposited a processed version of these pairs files in ftp://cbsuftp.tc.cornell.edu/danko/hub/MicroC_pairs_files/. You will also need the baits and prey files which can be found at the input files directory in this GitHub repository. To obtain the raw aggregated contacts 20kbX20kb matrix at 200bp resolution between enhancers and promoters within 25-150kb of genomic distance, as well as the vectors for the 1D signal around these enhancers and promoters at the DMSO control, use:

`bash MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh mESCs_DMSO_30_intra.mm10.nodups.pairs.gz dREG_based_promoters_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed dREG_based_TREs_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed 25000 150000 10000 50 outputPath`

And for flavopiridol (FLV) treated cells: 

`bash MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh mESCs_FLV_30_intra.mm10.nodups.pairs.gz dREG_based_promoters_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed dREG_based_TREs_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed 25000 150000 10000 50 outputPath`

And for Triptolide (TRP) treated cells: 

`bash MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh mESCs_TRP_30_intra.mm10.nodups.pairs.gz dREG_based_promoters_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed dREG_based_TREs_with_STARTseq_based_maxTSS_mm10_200bp_centered_on_maxTSS_chr_start_end_strand.bed 25000 150000 10000 50 outputPath`

Note: outputPath refers to the directory where work will be done. The acompaning python files, MicroC_Stranded_Aggregation_pipeline_get_bait_matrix.py, MicroC_Stranded_Aggregation_pipeline_get_aggregated_matrix.py and get_genome_wide_normalization_scores_by_search_window.py should be at the same directory as ContactCaller_microC.bsh.

After running MicroC_Stranded_Aggregation_pipeline_with_1D_signal.bsh for all treatment and control samples, you can visualize the change APAs for 10530 (promoters) baits and 27900 (enhancers) preys as follows:

`python Change_calculation_and_visualization.py Control_observed_contacts_APA_matrix.csv Treatment_observed_contacts_APA_matrix.csv Control_baits_1D_signal_vector Control_preys_1D_signal_vectorTreatment_baits_1D_signal_vector Treatment_preys_1D_signal_vector 10000 50 10530 27900 1D_normalized_change_APA.svg` 
