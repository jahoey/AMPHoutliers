# AMPHoutliers
code for identifying locus-environmental associations in *Amphiprion clarkii*
- **clownfishoutliers.R** script generates files for BayEnv2, reads in and concatinates BayEnv2 analysis files & then identifies potential candidates
- **clownfishRDA.R** script conducts redundancy analysis & identifies candidate loci associated with the environment
## ABC folder: files and scripts for performing ABC
A guide to performing ABC can be found [here](https://github.com/jahoey/pinskylab_methods/blob/master/genomics/Guide%20to%20performing%20ABC%20using%20SNP%20data.md).
- **fsc_SFSmanipulation_APCL1_0.R** manipulates and downsamples raw sfs files from fastsimcoal, results in **sub.sfs.all1.RData**
- **fsc_SFSmanipulation_APCL2_0.R** manipulates and downsamples raw sfs files from fastsimcoal, results in **sub.sfs.all2.RData**
- **fsc_SFSmanipulation_APCL2_1.R** manipulates and downsamples raw sfs files from fastsimcoal, results in **sub.sfs.all3.RData**
- **sim_joined_sfs_3pops.RData** R data file containing **sub.sfs.all1.RData**, **sub.sfs.all2.RData**, and **sub.sfs.all3.RData** column bound together
- **obs.3sfs.txt** observed SFSs for 3 populations of clownfish
- **readparams_APCL.R** this scripts creates **APCLfull_params.txt**
- **APCLfull_params.txt** parameters file necessary for ABC
- **migrationABC_600.est** and **migrationABC_600.tpl** are input files necessary for generating simulations with fastsimcoal
- **ABC join.R** this script prepares files for ABC and performs ABC analysis
