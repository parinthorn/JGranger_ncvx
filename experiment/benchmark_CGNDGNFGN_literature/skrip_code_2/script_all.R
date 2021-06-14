source("0_library_declare.r")

source("3_TWO_STAGE.r")

source("TWO_STAGE_ourformat.r")

filename <- file.choose()
file <- readRDS(filename)

write.csv(file,"A_true.csv")