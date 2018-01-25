setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Clownfish outlier analysis")

# Reading in names of polymorphic loci
locnames <- read.table("output.hicov2.snps.only.poly_locnames.txt")
locnames.first <- locnames[!duplicated(locnames[,1]),] # names of the loci I want to keep in the vcf

# Read in vcf but skip all the comments. This contains polymorphic and monomorphic loci.
vcf <- read.table("output.hicov2.snps.only.vcf", skip = 103549, sep = "\t")
dim(vcf) # 5729 x 34 (these are polymorphic and non-polymorphic snps & this number matches with "output.hicov2.snps.only.txt")

# Keeping only the first SNP at each contig
vcf_firstsnps <- vcf[!duplicated(vcf[,1]),] # these are the first snps at each contig
dim(vcf_firstsnps) # 1135 x 34

# Now remove non-polymorphic SNPs whose locus names are not in locnames.first
vcf_firstsnps_poly <- vcf_firstsnps[vcf_firstsnps$V1 %in% locnames.first$V1, ]
dim(vcf_firstsnps_poly)

# Are the locus names the same?
vcf_firstsnps_poly[111,1]
locnames.first[111,1]

write.table(vcf_firstsnps_poly, "vcf4str.vcf", sep="\t", col.names = FALSE, row.names = FALSE)
# Manually removed all the " and replaced with nothing

#### Now that the STRUCTURE-formatted snp file is made, can start from here ####
library(ade4)
library(adegenet)
library(devtools)
library("hierfstat")
library(pegas)
library(fields)

# Reading in SNP data file containing only the first SNP at each locus
# clowns <- read.structure("clownfish_Dec14.str",
#                          n.ind = 25, n.loc = 1135, col.lab = 1, col.pop = 2, row.marknames = 1, 
#                          onerowperind = FALSE)

clowns <- read.structure("clownfish_Dec15.str",
                         n.ind = 25, n.loc = 1131, col.lab = 1, col.pop = 2, row.marknames = 1, 
                         onerowperind = FALSE)

which(clowns@loc.n.all > 2) # which snps have more than 2 alelles? None
which(clowns@loc.n.all < 2) # 4 loci are mono - allelic: SNPs 284, 734, 740, 965 (clownfish_Dec14.str)

clowns.nonas <- scaleGen(clowns, center = TRUE, scale = FALSE) # Centering the data. There is no need to fill in missing data because there is none
hist(clowns.nonas)
sum(is.na(clowns.nonas))

# Removing mono-allelic loci (clownfish_Dec15.str)
# clowns.nonas <- clowns.nonas[, -c(567, 1466, 1477, 1926)]

# Keeping only one allele per locus
clowns.nonas.25.1131 <- clowns.nonas[, seq(1, ncol(clowns.nonas),by = 2)] # including every other column
# clowns.nonas.25.1131 <- clowns.nonas[, -seq(1, ncol(clowns.nonas),by = 2)] # every second allele at each locus
dim(clowns.nonas.25.1131)

#### Read in the environmental data ####
envi <- read.csv("Sample_data_w_env.csv", header = TRUE)
dim(envi)
envi.sub <- envi[, c("Sample", "Country", "Latitude", "sss_mean", "sst_mean", "sst_max", "sst_min")]
envi.sub

# Order the environmental data so that individuals are in the same order as the genetic matrix - could not for the life of me figure out how to do this more elegantly...
# rn <- c("J11H", "J13H", "J15H", "J17H", "J19H", "J3H",  "J5H",  "J9H",  "N1H",  "N2H",  "N3H",  "N4H",  "N5H",  "N6H",  "N7H",  "P1H",  "P10H", "P14H", "P16H", "P18H", "P20H", "P22H", "P24H", "P3H",  "P6H" )
rn <- c(4,  5,  6,  7,  8,  1,  2,  3, 19, 20, 21, 22, 23, 24, 25, 9, 12, 13, 14, 15, 16, 17, 18, 10, 11)
envi.ordered <- envi.sub[c(rn),]

rownames(clowns.nonas.25.1131)[22] # same order? let's test
envi.ordered[22,1]

# Scale the environmental data
envi.scaled <- scale(data.matrix(envi.ordered[,3:7]))

# Performing full RDA
envi.scaled <- as.data.frame(envi.scaled)

library(vegan)

clowns.rda <- rda(clowns.nonas.25.1131 ~ Latitude + sss_mean + sst_mean + sst_max + sst_min, envi.scaled, scale = FALSE)
clowns.rda
summary(clowns.rda)
RsquareAdj (clowns.rda)
R2a.adults.rda <- RsquareAdj(clowns.rda)$adj.r.squared
coef(clowns.rda)
plot(clowns.rda, scaling = 3)
plot(clowns.rda, choices = c(1,2), scaling=3)
plot(clowns.rda, choices = c(1,3), scaling=3)
plot(clowns.rda, choices = c(2,3), scaling=3)
anova(clowns.rda)
anova(clowns.rda, by="axis", step=1000)

# Here are the scores for each locus
spp.scr <- scores(clowns.rda, display = "species", scaling = 3, choices = c(1,2,3))

# Forester et al. (2016) used loci with scores +/- 3 SD from mean score for that axis to ID outlier loci. Only used first 3 axes.
mean.rda <- colMeans(spp.scr)
sd.rda1 <- sd(spp.scr[,"RDA1"])
sd.rda2 <- sd(spp.scr[,"RDA2"])
sd.rda3 <- sd(spp.scr[,"RDA3"])

rda1.hi <- mean.rda[1] + 3*sd.rda1
rda1.lo <- mean.rda[1] - 3*sd.rda1
which(spp.scr[,"RDA1"] > rda1.hi)
which(spp.scr[,"RDA1"] < rda1.lo)

rda2.hi <- mean.rda[2] + 3*sd.rda2
rda2.lo <- mean.rda[2] - 3*sd.rda2
which(spp.scr[,"RDA2"] > rda2.hi)
which(spp.scr[,"RDA2"] < rda2.lo)

rda3.hi <- mean.rda[3] + 3*sd.rda3
rda3.lo <- mean.rda[3] - 3*sd.rda3
which(spp.scr[,"RDA3"] > rda3.hi)
which(spp.scr[,"RDA3"] < rda3.lo)

rda.cans <- c((which(spp.scr[,"RDA1"] > rda1.hi)),
              (which(spp.scr[,"RDA1"] < rda1.lo)),
              (which(spp.scr[,"RDA2"] > rda2.hi)),
              (which(spp.scr[,"RDA2"] < rda2.lo)),
              (which(spp.scr[,"RDA3"] > rda3.hi)),
              (which(spp.scr[,"RDA3"] < rda3.lo)))

# Be careful about the difference between SNP number and index number (not always the same because if mono - allelic loci were removed)              
rda.cans <- unique(sort(rda.cans)) # [64 loci] 20   51   97  103  121  153  158  169  225  235  260  286  311  316  317  327  359  391  408  409  443
                                  # 500  516  521  532  552  553  560  587  595  607 631  679  682  714  725  727  749  810  811  826  829  830  853  869  884  886  891  908  913  922  
                                  # 925  927  933  960  964  993 1026 1061 1063 1074 1076 1083 1109

# This commented out section is for clownfish_Dec14
# rda.can.snpsno <- c(158, 169, 225, 287, 317, 444, 728, 928, 1067, 20, 51, 103, 153, 360, 409, 410, 533, 554, 608, 632, 889, 894, 916, 925, 
#                      997, 121, 235, 553, 561, 715, 813, 930, 936, 968, 1113, 328, 501, 522, 588, 829, 832, 856, 887, 911, 963, 1087, 392, 683, 
#                      97, 260, 312, 318, 517, 596, 680, 726, 752, 814, 833, 872, 1030, 1065, 1078, 1080)
# 
# rda.can.snpsno <- unique(sort(rda.can.snpsno)) # 20.0   51.0   97.0  103.0  121.0  153.0  158.0  169.0  225.0  235.0  260.0  287.0  312.0  317.0  
#       # 318.0  328.0  360.0  392.0  409.0  410.0  444.0  501.0 517.0  522.0  533.0  553.0  554.0  561.0  588.0  596.0  608.0  632.0  680.0  683.0  715.0  
#       # 726  728.0  752.0  813.0  814.0  829.0  832.0  833.0  856.0 872.0  887.0  889.0  894.0  911.0  916.0  925.0  928.0  930.0  936.0  963.0  968.0  
#       # 997.0 1030.0 1065.0 1067.0 1078.0 1080.0 1087.0 1113.0

# Okay, so these are outlier SNPs, but how strong are their correlations with the environmental variables?
# Linear regression between allele frequencies and environmental variables for the RDA outliers
lat20 <- lm(clowns.nonas.25.1131[,20] ~ envi.scaled[,"Latitude"]) # 2.348581e-09
sss20 <- lm(clowns.nonas.25.1131[,20] ~ envi.scaled[,"sss_mean"]) # 2.958472e-05
sst20 <- lm(clowns.nonas.25.1131[,20] ~ envi.scaled[,"sst_mean"]) # 2.876843e-12

lat51 <- lm(clowns.nonas.25.1131[,51] ~ envi.scaled[,"Latitude"]) # 7.812954e-09
sss51 <- lm(clowns.nonas.25.1131[,51] ~ envi.scaled[,"sss_mean"]) # 4.352194e-05
sst51 <- lm(clowns.nonas.25.1131[,51] ~ envi.scaled[,"sst_mean"]) # 2.774672e-11

lat97 <- lm(clowns.nonas.25.1131[,97] ~ envi.scaled[,"Latitude"]) 
sss97 <- lm(clowns.nonas.25.1131[,97] ~ envi.scaled[,"sss_mean"]) 
sst97 <- lm(clowns.nonas.25.1131[,97] ~ envi.scaled[,"sst_mean"]) 

lat103 <- lm(clowns.nonas.25.1131[,103] ~ envi.scaled[,"Latitude"]) 
sss103 <- lm(clowns.nonas.25.1131[,103] ~ envi.scaled[,"sss_mean"]) 
sst103 <- lm(clowns.nonas.25.1131[,103] ~ envi.scaled[,"sst_mean"]) 

lat121 <- lm(clowns.nonas.25.1131[,121] ~ envi.scaled[,"Latitude"]) 
sss121 <- lm(clowns.nonas.25.1131[,121] ~ envi.scaled[,"sss_mean"]) 
sst121 <- lm(clowns.nonas.25.1131[,121] ~ envi.scaled[,"sst_mean"]) 

lat153 <- lm(clowns.nonas.25.1131[,153] ~ envi.scaled[,"Latitude"]) # 1.426583e-07
sss153 <- lm(clowns.nonas.25.1131[,153] ~ envi.scaled[,"sss_mean"]) # 0.0001244421
sst153 <- lm(clowns.nonas.25.1131[,153] ~ envi.scaled[,"sst_mean"]) # 3.36641e-09

lat158 <- lm(clowns.nonas.25.1131[,158] ~ envi.scaled[,"Latitude"]) # 9.21022e-07
sss158 <- lm(clowns.nonas.25.1131[,158] ~ envi.scaled[,"sss_mean"]) # 0.0002679248
sst158 <- lm(clowns.nonas.25.1131[,158] ~ envi.scaled[,"sst_mean"]) # 5.212585e-08

lat169 <- lm(clowns.nonas.25.1131[,169] ~ envi.scaled[,"Latitude"]) 
sss169 <- lm(clowns.nonas.25.1131[,169] ~ envi.scaled[,"sss_mean"]) 
sst169 <- lm(clowns.nonas.25.1131[,169] ~ envi.scaled[,"sst_mean"]) 

lat225 <- lm(clowns.nonas.25.1131[,225] ~ envi.scaled[,"Latitude"]) # 0.0006465781
sss225 <- lm(clowns.nonas.25.1131[,225] ~ envi.scaled[,"sss_mean"]) 
sst225 <- lm(clowns.nonas.25.1131[,225] ~ envi.scaled[,"sst_mean"]) # 2.336883e-05

lat235 <- lm(clowns.nonas.25.1131[,235] ~ envi.scaled[,"Latitude"]) 
sss235 <- lm(clowns.nonas.25.1131[,235] ~ envi.scaled[,"sss_mean"]) 
sst235 <- lm(clowns.nonas.25.1131[,235] ~ envi.scaled[,"sst_mean"]) 

lat260 <- lm(clowns.nonas.25.1131[,260] ~ envi.scaled[,"Latitude"]) 
sss260 <- lm(clowns.nonas.25.1131[,260] ~ envi.scaled[,"sss_mean"]) 
sst260 <- lm(clowns.nonas.25.1131[,260] ~ envi.scaled[,"sst_mean"]) 

lat286 <- lm(clowns.nonas.25.1131[,286] ~ envi.scaled[,"Latitude"]) # 3.607174e-05
sss286 <- lm(clowns.nonas.25.1131[,286] ~ envi.scaled[,"sss_mean"]) 
sst286 <- lm(clowns.nonas.25.1131[,286] ~ envi.scaled[,"sst_mean"]) # 7.020161e-06

lat311 <- lm(clowns.nonas.25.1131[,311] ~ envi.scaled[,"Latitude"]) 
sss311 <- lm(clowns.nonas.25.1131[,311] ~ envi.scaled[,"sss_mean"]) 
sst311 <- lm(clowns.nonas.25.1131[,311] ~ envi.scaled[,"sst_mean"]) 

lat316 <- lm(clowns.nonas.25.1131[,316] ~ envi.scaled[,"Latitude"]) # 0.0007908307
sss316 <- lm(clowns.nonas.25.1131[,316] ~ envi.scaled[,"sss_mean"]) 
sst316 <- lm(clowns.nonas.25.1131[,316] ~ envi.scaled[,"sst_mean"]) # 0.0001062192

lat317 <- lm(clowns.nonas.25.1131[,317] ~ envi.scaled[,"Latitude"]) 
sss317 <- lm(clowns.nonas.25.1131[,317] ~ envi.scaled[,"sss_mean"]) 
sst317 <- lm(clowns.nonas.25.1131[,317] ~ envi.scaled[,"sst_mean"]) 

lat327 <- lm(clowns.nonas.25.1131[,327] ~ envi.scaled[,"Latitude"]) 
sss327 <- lm(clowns.nonas.25.1131[,327] ~ envi.scaled[,"sss_mean"]) 
sst327 <- lm(clowns.nonas.25.1131[,327] ~ envi.scaled[,"sst_mean"]) 

lat359 <- lm(clowns.nonas.25.1131[,359] ~ envi.scaled[,"Latitude"]) # 2.231215e-05
sss359 <- lm(clowns.nonas.25.1131[,359] ~ envi.scaled[,"sss_mean"]) 
sst359 <- lm(clowns.nonas.25.1131[,359] ~ envi.scaled[,"sst_mean"]) # 3.205457e-07

lat391 <- lm(clowns.nonas.25.1131[,391] ~ envi.scaled[,"Latitude"]) 
sss391 <- lm(clowns.nonas.25.1131[,391] ~ envi.scaled[,"sss_mean"]) 
sst391 <- lm(clowns.nonas.25.1131[,391] ~ envi.scaled[,"sst_mean"]) 

lat408 <- lm(clowns.nonas.25.1131[,408] ~ envi.scaled[,"Latitude"]) # 2.492698e-06
sss408 <- lm(clowns.nonas.25.1131[,408] ~ envi.scaled[,"sss_mean"]) # 1.103083e-05
sst408 <- lm(clowns.nonas.25.1131[,408] ~ envi.scaled[,"sst_mean"]) # 2.739911e-05

lat409 <- lm(clowns.nonas.25.1131[,409] ~ envi.scaled[,"Latitude"]) # 1.426583e-07
sss409 <- lm(clowns.nonas.25.1131[,409] ~ envi.scaled[,"sss_mean"]) # 0.0001244421
sst409 <- lm(clowns.nonas.25.1131[,409] ~ envi.scaled[,"sst_mean"]) # 3.36641e-09

lat443 <- lm(clowns.nonas.25.1131[,443] ~ envi.scaled[,"Latitude"]) # 3.607174e-05
sss443 <- lm(clowns.nonas.25.1131[,443] ~ envi.scaled[,"sss_mean"]) 
sst443 <- lm(clowns.nonas.25.1131[,443] ~ envi.scaled[,"sst_mean"]) # 7.020161e-06

lat500 <- lm(clowns.nonas.25.1131[,500] ~ envi.scaled[,"Latitude"]) 
sss500 <- lm(clowns.nonas.25.1131[,500] ~ envi.scaled[,"sss_mean"]) 
sst500 <- lm(clowns.nonas.25.1131[,500] ~ envi.scaled[,"sst_mean"]) 

lat516 <- lm(clowns.nonas.25.1131[,516] ~ envi.scaled[,"Latitude"]) 
sss516 <- lm(clowns.nonas.25.1131[,516] ~ envi.scaled[,"sss_mean"]) 
sst516 <- lm(clowns.nonas.25.1131[,516] ~ envi.scaled[,"sst_mean"]) 

lat521 <- lm(clowns.nonas.25.1131[,521] ~ envi.scaled[,"Latitude"]) 
sss521 <- lm(clowns.nonas.25.1131[,521] ~ envi.scaled[,"sss_mean"]) 
sst521 <- lm(clowns.nonas.25.1131[,521] ~ envi.scaled[,"sst_mean"]) 

lat532 <- lm(clowns.nonas.25.1131[,532] ~ envi.scaled[,"Latitude"]) # 3.185348e-06
sss532 <- lm(clowns.nonas.25.1131[,532] ~ envi.scaled[,"sss_mean"]) # 0.0001291054
sst532 <- lm(clowns.nonas.25.1131[,532] ~ envi.scaled[,"sst_mean"]) # 2.743583e-06

lat552 <- lm(clowns.nonas.25.1131[,552] ~ envi.scaled[,"Latitude"]) 
sss552 <- lm(clowns.nonas.25.1131[,552] ~ envi.scaled[,"sss_mean"]) 
sst552 <- lm(clowns.nonas.25.1131[,552] ~ envi.scaled[,"sst_mean"]) 

lat553 <- lm(clowns.nonas.25.1131[,553] ~ envi.scaled[,"Latitude"]) # 7.18745e-05
sss553 <- lm(clowns.nonas.25.1131[,553] ~ envi.scaled[,"sss_mean"]) # 0.0006358537
sst553 <- lm(clowns.nonas.25.1131[,553] ~ envi.scaled[,"sst_mean"]) # 9.488234e-05

lat560 <- lm(clowns.nonas.25.1131[,560] ~ envi.scaled[,"Latitude"]) 
sss560 <- lm(clowns.nonas.25.1131[,560] ~ envi.scaled[,"sss_mean"]) 
sst560 <- lm(clowns.nonas.25.1131[,560] ~ envi.scaled[,"sst_mean"]) 

lat587 <- lm(clowns.nonas.25.1131[,587] ~ envi.scaled[,"Latitude"]) 
sss587 <- lm(clowns.nonas.25.1131[,587] ~ envi.scaled[,"sss_mean"]) 
sst587 <- lm(clowns.nonas.25.1131[,587] ~ envi.scaled[,"sst_mean"]) 

lat595 <- lm(clowns.nonas.25.1131[,595] ~ envi.scaled[,"Latitude"]) 
sss595 <- lm(clowns.nonas.25.1131[,595] ~ envi.scaled[,"sss_mean"]) 
sst595 <- lm(clowns.nonas.25.1131[,595] ~ envi.scaled[,"sst_mean"]) 

lat607 <- lm(clowns.nonas.25.1131[,607] ~ envi.scaled[,"Latitude"]) 
sss607 <- lm(clowns.nonas.25.1131[,607] ~ envi.scaled[,"sss_mean"]) 
sst607 <- lm(clowns.nonas.25.1131[,607] ~ envi.scaled[,"sst_mean"]) # 0.0004916976

lat631 <- lm(clowns.nonas.25.1131[,631] ~ envi.scaled[,"Latitude"]) # 2.376022e-05
sss631 <- lm(clowns.nonas.25.1131[,631] ~ envi.scaled[,"sss_mean"]) 
sst631 <- lm(clowns.nonas.25.1131[,631] ~ envi.scaled[,"sst_mean"]) # 2.391972e-06

lat679 <- lm(clowns.nonas.25.1131[,679] ~ envi.scaled[,"Latitude"]) 
sss679 <- lm(clowns.nonas.25.1131[,679] ~ envi.scaled[,"sss_mean"]) 
sst679 <- lm(clowns.nonas.25.1131[,679] ~ envi.scaled[,"sst_mean"]) 

lat682 <- lm(clowns.nonas.25.1131[,682] ~ envi.scaled[,"Latitude"]) 
sss682 <- lm(clowns.nonas.25.1131[,682] ~ envi.scaled[,"sss_mean"]) 
sst682 <- lm(clowns.nonas.25.1131[,682] ~ envi.scaled[,"sst_mean"]) 

lat714 <- lm(clowns.nonas.25.1131[,714] ~ envi.scaled[,"Latitude"]) 
sss714 <- lm(clowns.nonas.25.1131[,714] ~ envi.scaled[,"sss_mean"]) 
sst714 <- lm(clowns.nonas.25.1131[,714] ~ envi.scaled[,"sst_mean"]) 

lat725 <- lm(clowns.nonas.25.1131[,725] ~ envi.scaled[,"Latitude"]) 
sss725 <- lm(clowns.nonas.25.1131[,725] ~ envi.scaled[,"sss_mean"]) 
sst725 <- lm(clowns.nonas.25.1131[,725] ~ envi.scaled[,"sst_mean"]) 

lat727 <- lm(clowns.nonas.25.1131[,727] ~ envi.scaled[,"Latitude"]) # 4.146615e-05
sss727 <- lm(clowns.nonas.25.1131[,727] ~ envi.scaled[,"sss_mean"]) 
sst727 <- lm(clowns.nonas.25.1131[,727] ~ envi.scaled[,"sst_mean"]) # 8.377551e-06

lat749 <- lm(clowns.nonas.25.1131[,749] ~ envi.scaled[,"Latitude"]) 
sss749 <- lm(clowns.nonas.25.1131[,749] ~ envi.scaled[,"sss_mean"]) 
sst749 <- lm(clowns.nonas.25.1131[,749] ~ envi.scaled[,"sst_mean"]) 

lat810 <- lm(clowns.nonas.25.1131[,810] ~ envi.scaled[,"Latitude"]) 
sss810 <- lm(clowns.nonas.25.1131[,810] ~ envi.scaled[,"sss_mean"]) 
sst810 <- lm(clowns.nonas.25.1131[,810] ~ envi.scaled[,"sst_mean"]) 

lat811 <- lm(clowns.nonas.25.1131[,811] ~ envi.scaled[,"Latitude"]) 
sss811 <- lm(clowns.nonas.25.1131[,811] ~ envi.scaled[,"sss_mean"]) 
sst811 <- lm(clowns.nonas.25.1131[,811] ~ envi.scaled[,"sst_mean"]) 

lat826 <- lm(clowns.nonas.25.1131[,826] ~ envi.scaled[,"Latitude"]) 
sss826 <- lm(clowns.nonas.25.1131[,826] ~ envi.scaled[,"sss_mean"]) 
sst826 <- lm(clowns.nonas.25.1131[,826] ~ envi.scaled[,"sst_mean"]) 

lat829 <- lm(clowns.nonas.25.1131[,829] ~ envi.scaled[,"Latitude"]) 
sss829 <- lm(clowns.nonas.25.1131[,829] ~ envi.scaled[,"sss_mean"]) 
sst829 <- lm(clowns.nonas.25.1131[,829] ~ envi.scaled[,"sst_mean"]) 

lat830 <- lm(clowns.nonas.25.1131[,830] ~ envi.scaled[,"Latitude"]) 
sss830 <- lm(clowns.nonas.25.1131[,830] ~ envi.scaled[,"sss_mean"]) 
sst830 <- lm(clowns.nonas.25.1131[,830] ~ envi.scaled[,"sst_mean"]) 

lat853 <- lm(clowns.nonas.25.1131[,853] ~ envi.scaled[,"Latitude"]) 
sss853 <- lm(clowns.nonas.25.1131[,853] ~ envi.scaled[,"sss_mean"]) 
sst853 <- lm(clowns.nonas.25.1131[,853] ~ envi.scaled[,"sst_mean"]) 

lat869 <- lm(clowns.nonas.25.1131[,869] ~ envi.scaled[,"Latitude"]) 
sss869 <- lm(clowns.nonas.25.1131[,869] ~ envi.scaled[,"sss_mean"]) 
sst869 <- lm(clowns.nonas.25.1131[,869] ~ envi.scaled[,"sst_mean"]) 

lat884 <- lm(clowns.nonas.25.1131[,884] ~ envi.scaled[,"Latitude"]) 
sss884 <- lm(clowns.nonas.25.1131[,884] ~ envi.scaled[,"sss_mean"]) 
sst884 <- lm(clowns.nonas.25.1131[,884] ~ envi.scaled[,"sst_mean"]) 

lat886 <- lm(clowns.nonas.25.1131[,886] ~ envi.scaled[,"Latitude"]) 
sss886 <- lm(clowns.nonas.25.1131[,886] ~ envi.scaled[,"sss_mean"]) 
sst886 <- lm(clowns.nonas.25.1131[,886] ~ envi.scaled[,"sst_mean"]) # 0.0005662562

lat891 <- lm(clowns.nonas.25.1131[,891] ~ envi.scaled[,"Latitude"]) 
sss891 <- lm(clowns.nonas.25.1131[,891] ~ envi.scaled[,"sss_mean"]) 
sst891 <- lm(clowns.nonas.25.1131[,891] ~ envi.scaled[,"sst_mean"]) # 0.0005662562

lat908 <- lm(clowns.nonas.25.1131[,908] ~ envi.scaled[,"Latitude"]) 
sss908 <- lm(clowns.nonas.25.1131[,908] ~ envi.scaled[,"sss_mean"]) 
sst908 <- lm(clowns.nonas.25.1131[,908] ~ envi.scaled[,"sst_mean"]) 

lat913 <- lm(clowns.nonas.25.1131[,913] ~ envi.scaled[,"Latitude"]) 
sss913 <- lm(clowns.nonas.25.1131[,913] ~ envi.scaled[,"sss_mean"]) 
sst913 <- lm(clowns.nonas.25.1131[,913] ~ envi.scaled[,"sst_mean"]) 

lat922 <- lm(clowns.nonas.25.1131[,922] ~ envi.scaled[,"Latitude"]) # 9.578872e-06
sss922 <- lm(clowns.nonas.25.1131[,922] ~ envi.scaled[,"sss_mean"]) # 0.0004997594
sst922 <- lm(clowns.nonas.25.1131[,922] ~ envi.scaled[,"sst_mean"]) # 2.859658e-06

lat925 <- lm(clowns.nonas.25.1131[,925] ~ envi.scaled[,"Latitude"]) # 5.7605e-05
sss925 <- lm(clowns.nonas.25.1131[,925] ~ envi.scaled[,"sss_mean"])
sst925 <- lm(clowns.nonas.25.1131[,925] ~ envi.scaled[,"sst_mean"]) # 3.202362e-06

lat927 <- lm(clowns.nonas.25.1131[,927] ~ envi.scaled[,"Latitude"]) 
sss927 <- lm(clowns.nonas.25.1131[,927] ~ envi.scaled[,"sss_mean"])
sst927 <- lm(clowns.nonas.25.1131[,927] ~ envi.scaled[,"sst_mean"]) 

lat933 <- lm(clowns.nonas.25.1131[,933] ~ envi.scaled[,"Latitude"]) 
sss933 <- lm(clowns.nonas.25.1131[,933] ~ envi.scaled[,"sss_mean"]) # 6.134071e-05
sst933 <- lm(clowns.nonas.25.1131[,933] ~ envi.scaled[,"sst_mean"]) 

lat960 <- lm(clowns.nonas.25.1131[,960] ~ envi.scaled[,"Latitude"]) 
sss960 <- lm(clowns.nonas.25.1131[,960] ~ envi.scaled[,"sss_mean"]) 
sst960 <- lm(clowns.nonas.25.1131[,960] ~ envi.scaled[,"sst_mean"]) 

lat964 <- lm(clowns.nonas.25.1131[,964] ~ envi.scaled[,"Latitude"]) 
sss964 <- lm(clowns.nonas.25.1131[,964] ~ envi.scaled[,"sss_mean"]) # 7.485169e-06
sst964 <- lm(clowns.nonas.25.1131[,964] ~ envi.scaled[,"sst_mean"]) 

lat993 <- lm(clowns.nonas.25.1131[,993] ~ envi.scaled[,"Latitude"]) # 0.0007781501
sss993 <- lm(clowns.nonas.25.1131[,993] ~ envi.scaled[,"sss_mean"]) 
sst993 <- lm(clowns.nonas.25.1131[,993] ~ envi.scaled[,"sst_mean"]) # 0.0003499796

lat1026 <- lm(clowns.nonas.25.1131[,1026] ~ envi.scaled[,"Latitude"]) 
sss1026 <- lm(clowns.nonas.25.1131[,1026] ~ envi.scaled[,"sss_mean"]) 
sst1026 <- lm(clowns.nonas.25.1131[,1026] ~ envi.scaled[,"sst_mean"]) 

lat1061 <- lm(clowns.nonas.25.1131[,1061] ~ envi.scaled[,"Latitude"]) 
sss1061 <- lm(clowns.nonas.25.1131[,1061] ~ envi.scaled[,"sss_mean"]) 
sst1061 <- lm(clowns.nonas.25.1131[,1061] ~ envi.scaled[,"sst_mean"]) 

lat1063 <- lm(clowns.nonas.25.1131[,1063] ~ envi.scaled[,"Latitude"]) # 1.819428e-06
sss1063 <- lm(clowns.nonas.25.1131[,1063] ~ envi.scaled[,"sss_mean"]) # 0.0003611427
sst1063 <- lm(clowns.nonas.25.1131[,1063] ~ envi.scaled[,"sst_mean"]) # 1.351008e-07

lat1074 <- lm(clowns.nonas.25.1131[,1074] ~ envi.scaled[,"Latitude"]) 
sss1074 <- lm(clowns.nonas.25.1131[,1074] ~ envi.scaled[,"sss_mean"]) 
sst1074 <- lm(clowns.nonas.25.1131[,1074] ~ envi.scaled[,"sst_mean"]) 

lat1076 <- lm(clowns.nonas.25.1131[,1076] ~ envi.scaled[,"Latitude"]) 
sss1076 <- lm(clowns.nonas.25.1131[,1076] ~ envi.scaled[,"sss_mean"]) 
sst1076 <- lm(clowns.nonas.25.1131[,1076] ~ envi.scaled[,"sst_mean"]) 

lat1083 <- lm(clowns.nonas.25.1131[,1083] ~ envi.scaled[,"Latitude"]) 
sss1083 <- lm(clowns.nonas.25.1131[,1083] ~ envi.scaled[,"sss_mean"]) 
sst1083 <- lm(clowns.nonas.25.1131[,1083] ~ envi.scaled[,"sst_mean"]) 

lat1109 <- lm(clowns.nonas.25.1131[,1109] ~ envi.scaled[,"Latitude"]) 
sss1109 <- lm(clowns.nonas.25.1131[,1109] ~ envi.scaled[,"sss_mean"]) 
sst1109 <- lm(clowns.nonas.25.1131[,1109] ~ envi.scaled[,"sst_mean"]) 

# Function to calculate p-values
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

pvals <- c(lmp(lat20),
           lmp(sss20),
           lmp(sst20),
           lmp(lat51),
           lmp(sss51),
           lmp(sst51),
           lmp(lat97),
           lmp(sss97),
           lmp(sst97),
           lmp(lat103),
           lmp(sss103),
           lmp(sst103),
           lmp(lat121),
           lmp(sss121),
           lmp(sst121),
           lmp(lat153),
           lmp(sss153),
           lmp(sst153),
           lmp(lat158),
           lmp(sss158),
           lmp(sst158),
           lmp(lat169),
           lmp(sss169),
           lmp(sst169),
           lmp(lat225),
           lmp(sss225),
           lmp(sst225),
           lmp(lat235),
           lmp(sss235),
           lmp(sst235),
           lmp(lat260),
           lmp(sss260),
           lmp(sst260),
           lmp(lat286),
           lmp(sss286),
           lmp(sst286),
           lmp(lat311),
           lmp(sss311),
           lmp(sst311),
           lmp(lat316),
           lmp(sss316),
           lmp(sst316),
           lmp(lat317),
           lmp(sss317),
           lmp(sst317),
           lmp(lat327),
           lmp(sss327),
           lmp(sst327),
           lmp(lat359),
           lmp(sss359),
           lmp(sst359),
           lmp(lat391),
           lmp(sss391),
           lmp(sst391),
           lmp(lat408),
           lmp(sss408),
           lmp(sst408),
           lmp(lat409),
           lmp(sss409),
           lmp(sst409),
           lmp(lat443),
           lmp(sss443),
           lmp(sst443),
           lmp(lat500),
           lmp(sss500),
           lmp(sst500),
           lmp(lat516),
           lmp(sss516),
           lmp(sst516),
           lmp(lat521),
           lmp(sss521),
           lmp(sst521),
           lmp(lat532),
           lmp(sss532),
           lmp(sst532),
           lmp(lat552),
           lmp(sss552),
           lmp(sst552),
           lmp(lat553),
           lmp(sss553),
           lmp(sst553),
           lmp(lat560),
           lmp(sss560),
           lmp(sst560),
           lmp(lat587),
           lmp(sss587),
           lmp(sst587),
           lmp(lat595),
           lmp(sss595),
           lmp(sst595),
           lmp(lat607),
           lmp(sss607),
           lmp(sst607),
           lmp(lat631),
           lmp(sss631),
           lmp(sst631),
           lmp(lat679),
           lmp(sss679),
           lmp(sst679),
           lmp(lat682),
           lmp(sss682),
           lmp(sst682),
           lmp(lat714),
           lmp(sss714),
           lmp(sst714),
           lmp(lat725),
           lmp(sss725),
           lmp(sst725),
           lmp(lat727),
           lmp(sss727),
           lmp(sst727),
           lmp(lat749),
           lmp(sss749),
           lmp(sst749),
           lmp(lat810),
           lmp(sss810),
           lmp(sst810),
           lmp(lat811),
           lmp(sss811),
           lmp(sst811),
           lmp(lat826),
           lmp(sss826),
           lmp(sst826),
           lmp(lat829),
           lmp(sss829),
           lmp(sst829),
           lmp(lat830),
           lmp(sss830),
           lmp(sst830),
           lmp(lat853),
           lmp(sss853),
           lmp(sst853),
           lmp(lat869),
           lmp(sss869),
           lmp(sst869),
           lmp(lat884),
           lmp(sss884),
           lmp(sst884),
           lmp(lat886),
           lmp(sss886),
           lmp(sst886),
           lmp(lat891),
           lmp(sss891),
           lmp(sst891),
           lmp(lat908),
           lmp(sss908),
           lmp(sst908),
           lmp(lat913),
           lmp(sss913),
           lmp(sst913),
           lmp(lat922),
           lmp(sss922),
           lmp(sst922),
           lmp(lat925),
           lmp(sss925),
           lmp(sst925),
           lmp(lat927),
           lmp(sss927),
           lmp(sst927),
           lmp(lat933),
           lmp(sss933),
           lmp(sst933),
           lmp(lat960),
           lmp(sss960),
           lmp(sst960),
           lmp(lat964),
           lmp(sss964),
           lmp(sst964),
           lmp(lat993),
           lmp(sss993),
           lmp(sst993),
           lmp(lat1026),
           lmp(sss1026),
           lmp(sst1026),
           lmp(lat1061),
           lmp(sss1061),
           lmp(sst1061),
           lmp(lat1063),
           lmp(sss1063),
           lmp(sst1063),
           lmp(lat1074),
           lmp(sss1074),
           lmp(sst1074),
           lmp(lat1076),
           lmp(sss1076),
           lmp(sst1076),
           lmp(lat1083),
           lmp(sss1083),
           lmp(sst1083),
           lmp(lat1109),
           lmp(sss1109),
           lmp(sst1109))

which(pvals < 0.001)

# Plotting the RDA plot with loci indicated by redundancy analysis to have strong locus-environmental associations
png(file="~/Documents/Graduate School/Rutgers/Clownfish outlier analysis/RDA plot.png", width=12, height=6, res=300, units="in")

par(
  mfrow = c(1, 2), 
  mar=c(5, 5, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg="white"
)

plot(clowns.rda, choices = c(1,2), scaling=3, type = "none")
points(clowns.rda, choices = c(1,2), display = "sites", col = "gray80", scaling=3, cex = 0.6)
points(clowns.rda, choices = c(1,2), display = "sp", col = "gray50", scaling=3, pch = 3, cex = 0.6)
points(clowns.rda, choices = c(1,2), display = "bp", scaling=3)
text(clowns.rda, choices = c(1,2), display = "bp", scaling=3, cex = 0.6)

spp.scr <- scores(clowns.rda, display = "species", scaling = 3, choices = c(1,2,3))
# site.scr <- scores(adults.rda, display = "sites", scaling = 3, choices = c(1,2,3,4,5))
# bp.scr <- scores(adults.rda, display = "bp", scaling = 3, choices = c(1,2,3,4,5))
points(spp.scr[20,1], spp.scr[20,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[51,1], spp.scr[51,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[153,1], spp.scr[153,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[158,1], spp.scr[158,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[225,1], spp.scr[225,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[286,1], spp.scr[286,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[316,1], spp.scr[316,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[359,1], spp.scr[359,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[408,1], spp.scr[408,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[409,1], spp.scr[409,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[443,1], spp.scr[443,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[532,1], spp.scr[532,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[553,1], spp.scr[553,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[607,1], spp.scr[607,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[631,1], spp.scr[631,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[727,1], spp.scr[727,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[886,1], spp.scr[886,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[891,1], spp.scr[891,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[922,1], spp.scr[922,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[925,1], spp.scr[925,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[933,1], spp.scr[933,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[964,1], spp.scr[964,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[993,1], spp.scr[993,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[1063,1], spp.scr[1063,2], col = "black", pch = 16, cex = 1.1)

plot(clowns.rda, choices = c(1,3), scaling=3, type = "none")
points(clowns.rda, choices = c(1,3), display = "sites", col = "gray80", scaling=3, cex = 0.6)
points(clowns.rda, choices = c(1,3), display = "sp", col = "gray50", scaling=3, pch = 3, cex = 0.6)
points(clowns.rda, choices = c(1,3), display = "bp", scaling=3)
text(clowns.rda, choices = c(1,3), display = "bp", scaling=3, cex = 0.6)
points(spp.scr[20,1], spp.scr[20,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[51,1], spp.scr[51,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[153,1], spp.scr[153,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[158,1], spp.scr[158,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[225,1], spp.scr[225,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[286,1], spp.scr[286,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[316,1], spp.scr[316,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[359,1], spp.scr[359,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[408,1], spp.scr[408,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[409,1], spp.scr[409,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[443,1], spp.scr[443,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[532,1], spp.scr[532,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[553,1], spp.scr[553,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[607,1], spp.scr[607,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[631,1], spp.scr[631,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[727,1], spp.scr[727,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[886,1], spp.scr[886,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[891,1], spp.scr[891,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[922,1], spp.scr[922,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[925,1], spp.scr[925,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[933,1], spp.scr[933,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[964,1], spp.scr[964,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[993,1], spp.scr[993,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[1063,1], spp.scr[1063,3], col = "black", pch = 16, cex = 1.1)

dev.off()

