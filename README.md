# admixture_load_scripts
Scripts for simulating the admixture models described in Kim et al. This is a brief overview. Please see the paper for specific details on what we did.

These shell scripts generate SLiM simulation scripts (*.slim files) which can then be used to recreate the simulated datasets in our paper. Because I ran these scripts on our cluster and wanted to avoid issues with random number generation for job arrays, I used a job ID plus a random integer as the random number seed for each simulation replicate. The simulations are run with the following command:

slim -seed ${JOB_ID} slimscript.slim

Most of these scripts can work with SLiM 2.x, but we recommend using SLiM 3.0 for full compatibility. Tree sequence recording allows us to accurately track ancestry and makes simulations much faster. In addition, all simulations are scaled to some factor which is defined in the simulation script files.

Msprime, pyslim, and pandas are required to run the Python scripts that parse the tree sequence output.

## 5Mb_simulations

This folder contains the 5Mb chunks that are being simulated. This outputs a CSV file which tracks population genetic statistics of interest over time.

## human_chr1_simulations

This folder contains simulation code for 100Mb of human chromosome 1. This outputs tree sequence files from which we parse out the local ancestry proportions in the genome.

The recombination map and exon definitions used for the simulations are found in sim_seq_info.txt.

## arabidopsis_selfing

This folder contains simulation code for 29.1 Mb of Arabidopsis chromosome 1. This outputs VCF files from which we parse out the local ancestry proportions in the genome. We did not use tree sequence recording because it was prohibitively slow due to higher recombination rates, i.e. there were more trees to record.

The recombination map and exon definitions used for the simulations are found in sim_seq_info.txt.