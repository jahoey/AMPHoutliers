setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Clownfish outlier analysis")

#### Creating a BayEnv formatted matrix to create the covariance matrix
# Read in allele frequency data split out by populations
snps.poly <- read.table("output.hicov2.snps.only.poly.txt") #only polymorphic snps
dim(snps.poly) #11436     3 (5718 snps)
str(snps.poly)
snps.poly["snp"] <- rep(1:5718, each = 2, length.out=11436) # adding a column for snp id

# Reading in loci names because I only want the first SNP at each contig
locnames <- read.table("output.hicov2.snps.only.poly_locnames.txt")
dim(locnames)
locnames.first <- locnames[!duplicated(locnames[,1]),] # this is the first snp at each contig, with snp number
dim(locnames.first)
rownames(locnames.first)
locnames.first["snp"] <- rownames(locnames.first)
dim(locnames.first)

# Remove any SNPs that aren't the first SNP on a contig
exclu_snps <- locnames[duplicated(locnames[,1]),] # these are duplicated snps on a contig
exclu_snps["snp"] <- rownames(exclu_snps)
exclu_snp_list <- exclu_snps[,"snp"]

snps.poly.first <- snps.poly[ ! snps.poly$snp %in% exclu_snp_list, ] # only taking the first snp at each contig and excluding all the duplicates from exclu_snp_list
dim(snps.poly.first)

snps.poly.first[567:568,] # check to make sure snp #'s match
locnames.first[284,]
# snps.poly.first[,"snp"] <- NULL # use this for the full matrix of 1131 snps
# write.table(snps.poly.first, "1131snpsfile.txt", sep="\t", col.names= FALSE, row.names = FALSE) # writes full matrix of allele counts for 1131 snps

# Write first SNP at each contig to a .txt file on my computer
# write.table (snps.poly.first, "1131loci_forhwe.txt", sep ="\t", col.names = FALSE, row.names = FALSE)

# Subsample a large number of SNPs to create the covariance matrix, and do this multiple times
snp.names <- snps.poly.first[duplicated(snps.poly.first[,"snp"]),4] # list of 1131 snp names
random1 <- sample(snp.names, 900, replace = FALSE) # randomly subsample 900 snp names
random.ord1 <- random[order(random1)] # this is a list of 900 randomly sampled snps
match(random.ord1,snp.names) # yes, elements of random.ord are contained within snp.names
random.allelecounts1 <- snps.poly.first[snps.poly.first$snp %in% random.ord1, ] # subsetting 1131 allele counts to only contain the snps in random.ord1
dim(random.allelecounts1)
random.allelecounts1[,"snp"] <- NULL
write.table(random.allelecounts1, "forbayenv_900covar1.txt", sep="\t", col.names = FALSE, row.names = FALSE)

# Another random draw of 900 snps
random2 <- sample(snp.names, 900, replace = FALSE) # randomly subsample 900 snp names
random.ord2 <- random2[order(random2)] # this is a list of 900 randomly sampled snps
match(random.ord2,snp.names) # yes, elements of random.ord are contained within snp.names
random.ord2 == random.ord1 # make sure it's different from the first random draw
random.allelecounts2 <- snps.poly.first[snps.poly.first$snp %in% random.ord2, ] # subsetting 1131 allele counts to only contain the snps in random.ord2
dim(random.allelecounts2)
random.allelecounts2[,"snp"] <- NULL
write.table(random.allelecounts2, "forbayenv_900covar2.txt", sep="\t", col.names = FALSE, row.names = FALSE)

#### Create the covariance matrix on Amphiprion and then average the covariance matricies ####

#### Averaging the MCMC output matricies so that I can use this averaged matrix as input for the next step of BayEnv2
# Reading in the matrix.out file from BayEnv2 using the fuction that Ryan wrote 

read_jho <- function(file, n_rows, n_matrices, n_comment, n_top_comment=0, ...){
  stopifnot(file.exists(file)) # check that file exists
  
  data_read <- list()
  for(i in 1:n_matrices){
    line_start <- (i-1)*(n_comment+n_rows) + n_comment + n_top_comment
    data_read[[i]] <- scan(file, skip=line_start, nlines=n_rows, ...)
    data_read[[i]] <- as.numeric(data_read[[i]])
    data_read[[i]] <- matrix(data_read[[i]], nrow=n_rows, byrow=TRUE)
  }
  
  data_array1 <- simplify2array(data_read)
  
  return(data_array1)
}

data_array1 <- read_jho(file="matrix1.out", n_rows=3, n_matrices=200, n_comment=2, n_top_comment=13, what="character")

# Take average of last 40 matrices
data_array1_mean <- apply(data_array1[,,160:200], c(1,2), mean)
data_array1_mean

# Reading in the second covariance matrix output
data_array2 <- read_jho(file="matrix2.out", n_rows=3, n_matrices=200, n_comment=2, n_top_comment=13, what="character")

# Take average of last 40 matrices
data_array2_mean <- apply(data_array2[,,160:200], c(1,2), mean)
data_array2_mean

# Visualizing the covariance matrices to compare them
library(corrplot)

covar1 <- data_array1_mean
covar1

cor1 <- cov2cor(covar1)
cor1

covar2 <- data_array2_mean
covar2

cor2 <- cov2cor(covar2)
cor2

# Plotting the covariance and correlation matricies
# Correlation matrix should reflect Fst between populations
par(mfrow=c(1,2))
corrplot(covar1, method = "color", main="Covariance Matrix1", is.corr = FALSE, cl.lim = c(-0.00009, 0.005))
corrplot(cor1, method = "color", main="Correlation Matrix1", is.corr = TRUE, cl.lim = c(-0.02, 1.000000000))

corrplot(covar2, method = "color", main="Covariance Matrix2", is.corr = FALSE) #cl.lim = c(-0.0000008, 0.005))
corrplot(cor2, method = "color", main="Correlation Matrix2", is.corr = TRUE, cl.lim = c(-0.03, 1.000000000))
# These plots look pretty similar, so I'll pick one and write it to my computer

write.table(data_array1_mean, "1131covarmatrix", sep="\t", col.names = FALSE, row.names = FALSE)

#### Read in the environmental data ####
envi <- read.csv("Sample_data_w_env.csv", header = TRUE)
dim(envi)
envi.sub <- envi[, c("Sample", "Country", "Latitude", "sss_mean", "sst_mean", "sst_max", "sst_min")]
envi.sub

# Aggregate the individuals into populations and standardize
by <- list(envi.sub$Country)
pop.avg <- aggregate(envi.sub[,3:7], by = by, FUN = mean)
stan.pop.avg <- scale(pop.avg[,2:6])
rownames(stan.pop.avg) <- c("Indonesia", "Japan", "Philippines")
t.stan.pop.avg <- t(stan.pop.avg) # the order is Indonesia, Japan, Philippines

# Exporting the envirofile to my computer
write.table (t.stan.pop.avg, "25stanenvirofile.amph.txt", sep ="\t", col.names = TRUE, row.names = TRUE) # manually redid the tabbing afterwards; order needs to match that of SNP frequencies (Japan, Indonesia, Philippines)

#### Covariance matrix created using 900 random loci, testing all 1131 loci for association with 5 environmental variables ####
#### Reading in BayEnv output files ####
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Clownfish outlier analysis/results_avg900randcovarmatrix_5stanenviro")

library(abind)

run1 <- read.table("bf1.txt")
run1$V1 <- NULL
run2 <- read.table("bf2.txt")
run2$V1 <- NULL
run3 <- read.table("bf3.txt")
run3$V1 <- NULL
run4 <- read.table("bf4.txt")
run4$V1 <- NULL
run5 <- read.table("bf5.txt")
run5$V1 <- NULL
run6 <- read.table("bf6.txt")
run6$V1 <- NULL
run7 <- read.table("bf7.txt")
run7$V1 <- NULL
run8 <- read.table("bf8.txt")
run8$V1 <- NULL
run9 <- read.table("bf9.txt")
run9$V1 <- NULL
run10 <- read.table("bf10.txt")
run10$V1 <- NULL

full_array <- abind(run1, run2, run3, run4, run5, run6, run7, run8, run9, run10, along=3)

full_array_mean <- apply(full_array, c(1,2), mean)
full_array_mean

which(full_array_mean[,1] > 3) # 20   51  151  153  158  164  169  286  359  392  408  409  443  449  469  478  504  532  553  568 591  631  727  886  891  913  922  925  933  993  999 1016 1050 1063 1081 1086 1108 1127
which(full_array_mean[,2] > 3) # 20   51  153  286  398  408  409  432  443  444  449  477  532  657  670  687  727  886  922  933  964 1063 1074 1127
which(full_array_mean[,3] > 3) # 20   51   58  151  153  158  164  181  187  192  205  223  225  247  270  286  322  331  359  392  408  409  443  469  504  532  553  568  607  631  727  886 891  897  911  913  922  925  928  993  999 1005 1063 1081 1108 1127
which(full_array_mean[,4] > 3) # 20   38   51   58  103  151  153  158  164  181  187  188  192  205  223  225  247  270  286  331  359  392  408  409  443  469  499  504  532  553  556  568 607  631  727  884  886  891  897  911  922  925  928  993  999 1005 1063 1081
which(full_array_mean[,5] > 3) # 20   51   58  151  153  158  164  181  187  192  205  223  225  247  270  286  322  331  359  392  408  409  443  469  504  532  553  568  607  631  727  886 891  897  911  913  922  925  928  993  999 1005 1063 1081 1108 1127

full_array_median <- apply(full_array, c(1,2), median)
full_array_median

which(full_array_median[,1] > 3) #[17]  20   51  153  158  286  359  408  409  443  532  553  727  891  913  922  993 1063
which(full_array_median[,2] > 3) #[8]  20   51  153  408  409  532  727 1063
which(full_array_median[,3] > 3) #[24]  20   51  153  158  225  286  359  392  408  409  443  532  553  607  631  727  886  891  913  922  925  928  993 1063
which(full_array_median[,4] > 3) #[24]  20   51  103  153  158  225  286  359  392  408  409  443  532  553  607  631  727  886  891  922  925  928  993 1063
which(full_array_median[,5] > 3) #[21]  20   51  153  158  286  359  408  409  443  532  553  607  631  727  886  891  913  922  925  993 1063

length(unique(c(which(full_array_median[,1] > 3), which(full_array_median[,2] > 3), which(full_array_median[,3] > 3), which(full_array_median[,4] > 3), which(full_array_median[,5] > 3)))) #[25] 20   51  103  153  158  225  286  359  392  408  409  443  532  553  607  631  727  886  891  913  922  925  928  993 1063

list<-(which(full_array_median[,1] > 3))
list <- append(list,which(full_array_median[,2] > 3))
list <- append(list,which(full_array_median[,3] > 3))
list <- append(list,which(full_array_median[,4] > 3))
list <- append(list,which(full_array_median[,5] > 3))

# Now I want the ID names of the contigs containing candidate SNPs. I will use locnames.first from the beginning of the script #
dim(locnames.first)
can_loci <- unique(sort(list))
can_names <- locnames.first[c(can_loci),]
candidates <- cbind(count(list), can_names)

write.table(candidates, "clownfish_candidateloci", sep="\t", row.names = FALSE)

# To examine allele frequencies by population
snps.poly.first[snps.poly.first$"snp" == "65",]