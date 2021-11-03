# Read functions
source(".\\Nf_modeling\\random_sampling_inE.R")
source(".\\Nf_modeling\\fit_wn_maha_model.R")

# Packages


# Occurrence data and species IDs
# native range, used to fit the models
sp.occ.nat <- read.csv("./Nf_modeling/occ_native_range.csv",header=T)
# invasive range, used to evaluate the models
sp.occ.inv <- read.csv("./Nf_modeling/occ_invasive_range.csv",header=T)
# species IDs
(sp.ids <- unique(sp.occ.nat[,1]))
# sample sizes
n.occ <- table(as.factor(sp.occ.nat[,1]))

### END ####
# Laura Jimenez
# November, 2021