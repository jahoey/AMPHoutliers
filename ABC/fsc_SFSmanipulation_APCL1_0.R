# Before beginning, I need to append 200 files with 500 SFS together in a certain order
# To be done on Amphiprion. This code reads in multiple .obs files that are in the same directory
setwd("~/25_MolEcoABC/pop1_0") # This needs to be the directory where the Amarel SFS files are

file_list <- list.files(".", recursive = TRUE, pattern = ".obs$") # should contain 200 file names

####################################################################################
# Make sure I'm in the right directory
setwd("~/25_MolEcoABC/pop1_0") # This needs to be the directory where the Amarel SFS files are
# Set up dimensions within each clump of simulations
# =============
# = Read File =
# =============

# ---- Data Dimensions Needed for Read ----
# See "Setup" Section Above
header_lines <- 0
bonus_top_lines <- 2 # these would be gibberish comments to skip in each matrix
good_lines <- 21 # these are the data lines
n_matrices <- 500

# ---- Read in each data file [200] 1 Matrix at a Time [500] ----
full_read <- list() # Create empty list to hold each simulation chunk as they're read in

# Loop and Read
for(i in 1:length(file_list)){
  data_read <- list() # Create empty list for each simulation chunk
  for(j in 1:n_matrices){
    
    # Determine Line of File to Start Reading for this Iteration
    lines_read_previously <- (j-1)*(bonus_top_lines+good_lines)
    line_start <- header_lines + lines_read_previously + (bonus_top_lines)
    
    # Read in Desired Portion of File
    # The 'sep' argument should refer to what separates values in the numeric part 
    # See the `write.table(make_nm(), ...` piece of code under 'Add to File' section above
    data_read[[j]] <- scan(file_list[i], skip=line_start, nlines=good_lines, what = as.list(c("character", rep("integer",17))))[-1] # [-1] gets rid of row names
    data_read[[j]] <- lapply(data_read[[j]], as.integer) # if problems try as.numeric
    
    # Format to Matrix
    # The value returned from `scan()` is a vector
    # The number of columns is inferred; can be set explicitly for clarity, if known
    data_read[[j]] <- matrix(unlist(data_read[[j]]), nrow=good_lines)
  }
  full_read <- append(full_read, data_read) # should be 100000 x 357 (pop1_0)
}

# Format to array
data_array <- simplify2array(full_read) #dim should be 21 x 17 x 100000

# Save the data array so that I can load it quickly
save(data_array, file = "~/25_MolEcoABC/pop1_0/data_array.RData")
load(file = "~/25_MolEcoABC/pop1_0/data_array.RData")

# Creates a list of how many loci each simulatlion has - will be a lot because don't need to do the MAF step
loci <- list()

for(i in 1:length(data_array[1,1,])){
  loci[i] <- sum(data_array[,,i])
}

#### Downsample all the simulations ####
setwd("~/25_MolEcoABC/pop1_0")
load(file = "~/25_MolEcoABC/pop1_0/data_array.RData")
num_loci <- 600 # set number of loci to keep for each simulation
n_matrices <- 100000
sub.sfs.all1 <- data.frame()

for(i in 1:n_matrices){
  index.one <- which(data_array[,,i] == 1) # indices of elements exactly equal to 1
  index.greaterthanone <- which(data_array[,,i] > 1) # indices of elements greater than 1
  adds <- vector()
  for(j in 1:length(index.greaterthanone)){
    adds <- append(adds, (rep(index.greaterthanone[j], (data_array[,,i][index.greaterthanone][j])-1)))
  }
  full.nonzeros <- c(index.one, index.greaterthanone, adds)
  sub.sfs <- data.frame()
  sub.sfs <- as.data.frame(table(sample(full.nonzeros,num_loci)))
  
  array.index <- data.frame(Var1 = 1:length(as.list(data_array[,,i]))) # Creating a dataframe with the number of indices equal to the number of entries in a 21 x 17 matrix (357 entries)
  
  array.counts <- merge(array.index,sub.sfs, by = "Var1", all.x = TRUE) # Merge array.index with the randomly sampled index.counts so that I have a dataframe of how many counts at each index
  array.counts[is.na(array.counts)] <- 0 # Replace all the NA's with zeros
  
  sub.dat <- matrix(array.counts$Freq, nrow = 1, ncol = 357) # Format the list of counts into a matrix (nrow = 21, ncol = 17) or a single row (nrow = 1, ncol = 357); matrices are filled column-wise
  
  sub.sfs.all1 <- rbind(sub.sfs.all1, sub.dat)
}

# Save downsampled simulations as a R datatype
save(sub.sfs.all1, file = "~/25_MolEcoABC/pop1_0/sub.sfs.all1.RData")



#### Plot ####
# Convert dataframe to matrix
sfs <- data.matrix(sub.dat)

#### Plot matrix ####
min <- min(sfs, na.rm = TRUE)
max <- max(sfs,na.rm = TRUE)

rbPal <- colorRampPalette(c('blue', 'red'))
ColorLevels <- seq(min, max, length=length(rbPal(5)))

# Set layout.  We are going to include a colorbar next to plot.
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1),
       heights=c(1,1))

#plotting margins.  These seem to work well for me.
par(mar = c(5,7,2.5,1), font = 2)

# Plot it up!
image(1:ncol(sfs), 1:nrow(sfs), t(sfs),
      col=rbPal(5), xlab="Pop_1", ylab="Pop_2",
      axes=FALSE, zlim = c(min, max),
      main= NA, xlim = c(0.5, 158), ylim = c(0.5, 110))
abline(0,1)

# Now annotate the plot
box()
axis(side = 1, at=seq(1,158,1), labels=0:157,
     cex.axis=1.0)
axis(side = 2, at=seq(1,110,1), labels=0:109, las= 1,
     cex.axis=1)

# Add colorbar to second plot region
par(mar = c(3,2.5,2.5,2))
image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=rbPal(5),xlab="",ylab="",xaxt="n", las = 1)

#### Number of loci sampled ####
pop0 <- colSums(sfs) #need sfs as a dataframe
hist(pop0, breaks = 135)
plot(pop0)
lines(pop0)

pop1 <- rowSums(sfs)
hist(pop0, breaks = 97)
plot(pop1)
lines(pop1)

###########################################################
#### Often times the SNP counts are really hard to see because they're all concentrated in the bottom left corner. Here is the same plotting code to zoom in ####
#### Plot matrix ####
min <- min(sfs, na.rm = TRUE)
max <- max(sfs,na.rm = TRUE)

rbPal <- colorRampPalette(c('blue', 'red'))
ColorLevels <- seq(min, max, length=length(rbPal(50)))

# Set layout.  We are going to include a colorbar next to plot.
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1),
       heights=c(1,1))

#plotting margins.  These seem to work well for me.
par(mar = c(5,7,2.5,1), font = 2)

# Plot it up!
image(1:ncol(sfs), 1:nrow(sfs), t(sfs),
      col=rbPal(50), xlab="Pop_1", ylab="Pop_2",
      axes=FALSE, zlim = c(min, max),
      main= NA, xlim = c(0.5, 25), ylim = c(0.5, 25))

# Now annotate the plot
box()
axis(side = 1, at=seq(1,25,1), labels=0:24,
     cex.axis=1.0)
axis(side = 2, at=seq(1,25,1), labels=0:24, las= 1,
     cex.axis=1)

# Add colorbar to second plot region
par(mar = c(3,2.5,2.5,2))
image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=rbPal(50),xlab="",ylab="",xaxt="n", las = 1)

