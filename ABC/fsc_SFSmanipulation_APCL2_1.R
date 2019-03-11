# Before beginning, I need to append 200 files with 500 SFS together in a certain order
# To be done on Amphiprion. This code reads in multiple .obs files that are in the same directory
setwd("~/25_MolEcoABC/pop2_1") # This needs to be the directory where the Amarel SFS files are

file_list <- list.files(".", recursive = TRUE, pattern = ".obs$") # should contain 200 file names

####################################################################################
# Make sure I'm in the right directory
setwd("~/25_MolEcoABC/pop2_1") # This needs to be the directory where the Amarel SFS files are
# Set up dimensions within each clump of simulations
# =============
# = Read File =
# =============

# ---- Data Dimensions Needed for Read ----
# See "Setup" Section Above
header_lines <- 0
bonus_top_lines <- 2 # these would be gibberish comments to skip in each matrix
good_lines <- 15 # these are the data lines
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
    data_read[[j]] <- scan(file_list[i], skip=line_start, nlines=good_lines, what = as.list(c("character", rep("integer",21))))[-1] # [-1] gets rid of row names
    data_read[[j]] <- lapply(data_read[[j]], as.integer) # if problems try as.numeric
    
    # Format to Matrix
    # The value returned from `scan()` is a vector
    # The number of columns is inferred; can be set explicitly for clarity, if known
    data_read[[j]] <- matrix(unlist(data_read[[j]]), nrow=good_lines)
  }
  full_read <- append(full_read, data_read) # should be 100000 x 315 (pop2_1)
}

# Format to array
data_array <- simplify2array(full_read) #dim should be 15 x 21 x 100000

# Save the data array so that I can load it quickly
save(data_array, file = "~/25_MolEcoABC/pop2_1/data_array.RData")
load(file = "~/25_MolEcoABC/pop2_1/data_array.RData")

#### Downsample all the simulations ####
setwd("~/25_MolEcoABC/pop2_1")
load(file = "~/25_MolEcoABC/pop2_1/data_array.RData")
num_loci <- 600 # set number of loci to keep for each simulation
n_matrices <- 100000
sub.sfs.all3 <- data.frame()

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
  
  array.index <- data.frame(Var1 = 1:length(as.list(data_array[,,i]))) # Creating a dataframe with the number of indices equal to the number of entries in a 15 x 21 matrix (315 entries)
  
  array.counts <- merge(array.index,sub.sfs, by = "Var1", all.x = TRUE) # Merge array.index with the randomly sampled index.counts so that I have a dataframe of how many counts at each index
  array.counts[is.na(array.counts)] <- 0 # Replace all the NA's with zeros
  
  sub.dat <- matrix(array.counts$Freq, nrow = 1, ncol = 315) # Format the list of counts into a matrix (nrow = 15, ncol = 21) or a single row (nrow = 1, ncol = 315); matrices are filled column-wise
  
  sub.sfs.all3 <- rbind(sub.sfs.all3, sub.dat)
}

# Save downsampled simulations as a R datatype
save(sub.sfs.all3, file = "~/25_MolEcoABC/pop2_1/sub.sfs.all3.RData")