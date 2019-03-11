# This meant to be run on Amphiprion. Should have already moved all the .params files from slurm-out to Amphiprion

# The extra columns should have already also been removed. See cut_cols.sh

# Two ways to do this
setwd("~/25_MolEcoABC/params")

# Way 1
paramsfile_list <- list.files(".", pattern = "new") # should contain 200 file names
tables <- lapply(paramsfile_list, read.table, header = TRUE)
full_params <- do.call(rbind,tables) # data.frame

write.table(full_params, file = "APCLfull_params.txt", sep = "\t", row.names = FALSE)