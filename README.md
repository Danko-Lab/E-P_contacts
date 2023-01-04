# Code for Barshad et al. 2023
This repository includes the main pieces of code used to support the conclusion in the paper. We also provide input files used in the paper when possible and explain how to obtain the ones that were too big to be provided here.

In all of the provided code, the source for Micro-C contact information is from prosseced pairs files generated via the distiller-nf algorithm. Most of these processed files can be found in the GEO dataset for this project: GSE206131. Other are available at ftp://cbsuftp.tc.cornell.edu/danko/hub/. A simple way of obtaining these files is to run distiller-nf with: `parsing_options: '--add-columns mapq' drop_readid: True` and then run:
