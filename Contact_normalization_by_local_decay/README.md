This directory includes the required bash and python files to normalize contacts between given sets of bait and preys in the genome. 
For each bait-prey pair, this code will calculate the observed bait-prey contacts, based on the contacts falling within a defined window
around both bait and prey, as well as the expected contacts, given a LOWESS fit to the contact decay over genomic distance around the prey.
This LOWESS fit is further corrected to zero inflation by accounting to the distribution of zero prey-associated contacts at different
genomic distances around the prey. The results file will include the chromosome and bait and prey coordinates with the observed and the expected
values along with a Fisher's exact test p-value associated with observed contacts deviation from the expected value.

Usage:

    bash ContactCaller_microC.bsh bait_file prey_file contacts_file results_directory CPU DIST CAP

Required files and arguments:

bait_file - File of "bait" peaks to analyze (format: Chromosome <\t> bait-center-position)

prey_file - File of "prey" TSSs to analyze (format: Chromosome <\t> prey-center-position)

contacts_file - contact (pairs) file output from distiler-nf without pair name and with mapq column, which is filtered to keep only cis contacts with both sides having a mapq >30 using: 

    zcat perfix.pairs.gz | awk 'BEGIN {OFS = "\t"} ; {if ($1 == "." && $2 == $4 && $9 >= 30 && $10 >= 30) {print $2, $3, $4, $5, $6, $7, $8, $9, $10}}' > perfix.nodups_30_intra.pairs
results_file - Path to output file

CPU - [optional, default=30] indicating how many CPU cores can be used.

DIST - [optional, default=1000000] Half size of "prey" search window around "bait" position

CAP - [optional, default=2000] Half size of "bait" and "prey" contact capture window

Other requirements:

- The python file should be placed in the directory from which the bash script is being executed.

