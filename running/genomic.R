# code for running per-year island attribution

library(RColorBrewer)
library(stringr)
library(dplyr)

# load in our dataset (relative to current working directory)
db_file <- "20150603.csv"
db <- read.csv(db_file)

# run the number of cases we want and so on

# options are:
# 1. Max number of alleles to impute
# 2. Source list
# 3. Output folder (defaults to current directory + "out_[imputed]_[sources]")

source_file <- "4_source_imputed"

sources <- read.csv(paste0(source_file, ".csv"), colClasses="character")
source_map <- as.numeric(sources$Number)
names(source_map) <- sources$DataSource

source_label_map <- unique(sources %>% select(Number, Label))
source_labels <- str_replace(source_label_map$Label, "\\\\n", "\n")
names(source_labels) <- source_label_map$Number

# input parameters
alleles_to_impute <- max(c(suppressWarnings(as.numeric(sources$Imputed)), 0), na.rm=T)

# model fitting control
seeds     <- c(5) #,7,11,13,17)
num_iters <- 500
thinning  <- 50

human   <- "Human"
mlst_cols  <- c("ST", "ASP", "GLN", "GLT", "GLY", "PGM", "TKT", "UNC")

# Setup data
db <- db %>% filter(Imputed <= alleles_to_impute)

humans <- db %>% filter(Source == human)

animals <- db %>% filter(Source %in% sources$DataSource)
animals <- animals %>% mutate(Source = source_map[as.character(Source)])
animals <- animals %>% select(one_of(c(mlst_cols, "Source")))

datasets <- list()
datasets[["all"]] <- humans %>% filter(Year >= 2005 & Year <= 2014)

# make the temp directory
temp_dir <- file.path("temp", source_file)
dir.create(temp_dir, showWarnings=F, recursive=T)

for (year in 1:length(datasets)) {

  # make the year directory
  year_dir <- file.path(temp_dir, names(datasets)[year])
  dir.create(year_dir, showWarnings=F)

  # create human dataset (we use the same animal dataset each time)
  h_year <- datasets[[year]] %>% select(one_of(mlst_cols))
  h_year <- h_year %>% mutate(Source = 0)

  data <- rbind(h_year, animals)

  # save data file
  island_in <- "input.txt"
  write.table(data, file=file.path(year_dir, island_in), sep="\t", row.names=F)

  # copy our isource executable there (as isource is braindead about output shit)
  file.copy("../isource", year_dir, overwrite=TRUE)
  current_dir <- getwd()

  for (seed in seeds) {
    island_seed_out <- paste0("output_", seed, ".txt")
    command    <- paste0("./isource ", island_in, " ", island_seed_out, " ", num_iters, " ", thinning, " 1 -", seed)

    cat("---------------------------------------\n")
    cat("running island model with seed", seed, "\n")
    cat(command, "\n")
    cat("---------------------------------------\n")

    setwd(year_dir)
    system(command)
    setwd(current_dir)
  }

  # TODO: analyse the output
}
