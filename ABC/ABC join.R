# This code is to be run on Amphiprion. 

# Read in each of the observed jointSFSs
setwd('~/25_MolEcoABC')
sfs1 <- read.table('output.hicov2.snps.only.arl_jointMAFpop2_0.obs', skip = 1) # simulated pop1_0
sfs2 <- read.table('output.hicov2.snps.only.arl_jointMAFpop1_0.obs', skip = 1) # simulated pop2_0
sfs3 <- read.table('output.hicov2.snps.only.arl_jointMAFpop2_1.obs', skip = 1) # simulated pop2_1

# Round down so that subsampling is easier (1.5 goes to 1)
sfs1 <- floor(sfs1)
sfs2 <- floor(sfs2)
sfs3 <- floor(t(sfs3))

# Downsample each observed jointSFSs
# Pop1_0
num_loci <- 600 # set number of loci to keep for each simulation
n_matrices <- 1
sub.sfs.all1 <- data.frame()

index.one <- which(sfs1 == 1) # indices of elements exactly equal to 1
index.greaterthanone <- which(sfs1 > 1) # indices of elements greater than 1
adds <- vector()
for(j in 1:length(index.greaterthanone)){
  adds <- append(adds, (rep(index.greaterthanone[j], unlist(sfs1)[index.greaterthanone][j]-1)))
}
full.nonzeros <- c(index.one, index.greaterthanone, adds)
sub.sfs <- data.frame()
sub.sfs <- as.data.frame(table(sample(full.nonzeros,num_loci)))

array.index <- data.frame(Var1 = 1:length(unlist(sfs1))) # Creating a dataframe with the number of indices equal to the number of entries in a 21 x 17 matrix (357 entries)

array.counts <- merge(array.index,sub.sfs, by = "Var1", all.x = TRUE) # Merge array.index with the randomly sampled index.counts so that I have a dataframe of how many counts at each index
array.counts[is.na(array.counts)] <- 0 # Replace all the NA's with zeros

sub.dat <- matrix(array.counts$Freq, nrow = 1, ncol = 357) # Format the list of counts into a matrix (nrow = 21, ncol = 17) or a single row (nrow = 1, ncol = 357); matrices are filled column-wise

sub.sfs.all1 <- rbind(sub.sfs.all1, sub.dat)

# Save downsampled simulations as a R datatype
save(sub.sfs.all1, file = "~/25_MolEcoABC/obssfspop1_0_sub.RData")

# Pop2_0
num_loci <- 600 # set number of loci to keep for each simulation
n_matrices <- 1
sub.sfs.all2 <- data.frame()

index.one <- which(sfs2 == 1) # indices of elements exactly equal to 1
index.greaterthanone <- which(sfs2 > 1) # indices of elements greater than 1
adds <- vector()
for(j in 1:length(index.greaterthanone)){
  adds <- append(adds, (rep(index.greaterthanone[j], unlist(sfs2)[index.greaterthanone][j]-1)))
}
full.nonzeros <- c(index.one, index.greaterthanone, adds)
sub.sfs <- data.frame()
sub.sfs <- as.data.frame(table(sample(full.nonzeros,num_loci)))

array.index <- data.frame(Var1 = 1:length(unlist(sfs2))) # Creating a dataframe with the number of indices equal to the number of entries in a 15 x 17 matrix (255 entries)

array.counts <- merge(array.index,sub.sfs, by = "Var1", all.x = TRUE) # Merge array.index with the randomly sampled index.counts so that I have a dataframe of how many counts at each index
array.counts[is.na(array.counts)] <- 0 # Replace all the NA's with zeros

sub.dat <- matrix(array.counts$Freq, nrow = 1, ncol = 255) # Format the list of counts into a matrix (nrow = 15, ncol = 17) or a single row (nrow = 1, ncol = 255); matrices are filled column-wise

sub.sfs.all2 <- rbind(sub.sfs.all2, sub.dat)

# Save downsampled simulations as a R datatype
save(sub.sfs.all2, file = "~/25_MolEcoABC/obssfspop2_0_sub.RData")

# Pop2_1
num_loci <- 600 # set number of loci to keep for each simulation
n_matrices <- 1
sub.sfs.all3 <- data.frame()

index.one <- which(sfs3 == 1) # indices of elements exactly equal to 1
index.greaterthanone <- which(sfs3 > 1) # indices of elements greater than 1
adds <- vector()
for(j in 1:length(index.greaterthanone)){
  adds <- append(adds, (rep(index.greaterthanone[j], unlist(sfs3)[index.greaterthanone][j]-1)))
}
full.nonzeros <- c(index.one, index.greaterthanone, adds)
sub.sfs <- data.frame()
sub.sfs <- as.data.frame(table(sample(full.nonzeros,num_loci)))

array.index <- data.frame(Var1 = 1:length(unlist(sfs3))) # Creating a dataframe with the number of indices equal to the number of entries in a 15 x 21 matrix (315 entries)

array.counts <- merge(array.index,sub.sfs, by = "Var1", all.x = TRUE) # Merge array.index with the randomly sampled index.counts so that I have a dataframe of how many counts at each index
array.counts[is.na(array.counts)] <- 0 # Replace all the NA's with zeros

sub.dat <- matrix(array.counts$Freq, nrow = 1, ncol = 315) # Format the list of counts into a matrix (nrow = 15, ncol = 21) or a single row (nrow = 1, ncol = 315); matrices are filled column-wise

sub.sfs.all3 <- rbind(sub.sfs.all3, sub.dat)

# Save downsampled simulations as a R datatype
save(sub.sfs.all3, file = "~/25_MolEcoABC/obssfspop2_1_sub.RData")

#### Need to create a single vector by using cbind to join each individual SFS ####
setwd('~/25_MolEcoABC')
load(file = "~/25_MolEcoABC/obssfspop1_0_sub.RData") # sub.sfs.all1
load(file = "~/25_MolEcoABC/obssfspop2_0_sub.RData") # sub.sfs.all2
load(file = "~/25_MolEcoABC/obssfspop2_1_sub.RData") # sub.sfs.all3

joined_sfs_3pops <- cbind(sub.sfs.all1,sub.sfs.all2,sub.sfs.all3)
write.table(joined_sfs_3pops, 'obs.3sfs.txt')

#### Trim simulations so the number of simulations match the observed - already done during 'downsampling' bit (fsc_SFSmanipulation_APCL1_0.R, fsc_SFSmanipulation_APCL2_0.R, fsc_SFSmanipulation_APCL2_1.R). Next, read in downsampled SFSs, then cbind into a single vector ####
# On my computer
setwd('~/Documents/Graduate School/Rutgers/Clownfish outlier analysis/ABC/')
load(file = "sub.sfs.all1.RData") # sub.sfs.all1
load(file = "sub.sfs.all2.RData") # sub.sfs.all2
load(file = "sub.sfs.all3.RData") # sub.sfs.all3

# Match the number of simulations
# sub.sfs.all <- sub.sfs.all[1:730,]
# sub.sfs.all2 <- sub.sfs.all2[1:730,]
# sub.sfs.all3 <- sub.sfs.all3[1:730,]

sim_joined_sfs_3pops <- cbind(sub.sfs.all1, sub.sfs.all2, sub.sfs.all3)

# write.table(sub.sfs.all, "sub.sfs.allpop1_0.txt")
# write.table(sub.sfs.all2, "sub.sfs.allpop2_0.txt")
# write.table(sub.sfs.all3, "sub.sfs.allpop2_1.txt")

write.table(sim_joined_sfs_3pops, "sim.3sfs.txt")
save(sim_joined_sfs_3pops, file = "~/Documents/Graduate School/Rutgers/Clownfish outlier analysis/ABC/sim_joined_sfs_3pops.RData")

#### ABC ####
library(abc)
library(KScorrect)

# setwd('~/25_MolEcoABC')
setwd('~/Documents/Graduate School/Rutgers/Clownfish outlier analysis/ABC/')
params <- read.table('APCLfull_params.txt', header = TRUE)
# params <- params[1:730,]
obs <- read.table('obs.3sfs.txt')
load(file = "sim_joined_sfs_3pops.RData") # sim_joined_sfs_3pops

calcs <- abc(target = obs, param = params, sumstat = sim_joined_sfs_3pops, tol = 0.05, method = 'rejection')
summary(calcs)

par(mfrow = c(2,2))
plot(density(log10(calcs$unadj.values[,1])), col='red', xlab = 'log10(POPONE)', main = '')
plot(density(log10(calcs$unadj.values[,2])), col='red', xlab = 'log10(POPTWO)', main = '')
plot(density(log10(calcs$unadj.values[,3])), col='red', xlab = 'log10(POPTHREE)', main = '')
plot(density(calcs$unadj.values[,4]), col='red', xlab = 'DISP', main = '')

