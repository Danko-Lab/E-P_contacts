# Code for Barshad et al. 2023
This repository includes the main pieces of code used to support the conclusion in the paper. We also provide input files used in the paper when possible and explain how to obtain the ones that were too big to be provided here. The three main types of analyses are: Contact_normalization_by_local_decay, APA_and_inter-sample_APA and EP_contacts_compared_to_local_background. In each of the above repositories you can find a detailed explanation about the different scripts, what they do and requred\optional parameters and input files. Here, we provide examples for code execution with data analysed in this paper.   

### Pairs files
In all of the provided code, the source for Micro-C contact information is from prosseced pairs files generated via the distiller-nf algorithm. Most of these processed files can be found in the GEO dataset for this project: GSE206131. Other are available at ftp://cbsuftp.tc.cornell.edu/danko/hub/MicroC_pairs_files/. A simple way of obtaining these files is to run distiller-nf with: `parsing_options: '--add-columns mapq' drop_readid: True`. You can find examples of raw Micro-C data processing in the "Micro-C_basic_processing" directory that containes YML configuration files used for distiller-nf. After obtaining pairs file for each replicate in your data, run:
`zcat perfix.rep1.pairs.gz perfix.rep2.pairs.gz perfix.rep3.pairs.gz ... | awk 'BEGIN {OFS = "\t"} ; {if ($1 == "." && $2 == $4 && $9 >= 30 && $10 >= 30) {print $2, $3, $4, $5, $6, $7, $8, $9, $10}}' > perfix.nodups_30_intra.pairs`

### Contact_normalization_by_local_decay

