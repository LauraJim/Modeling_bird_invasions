# agaper
# Error in `$<-.data.frame`(`*tmp*`, "foldID", value = c(2L, 2L, 3L, 2L,  )): 
#   replacement has 10019 rows, data has 10004
# 
# > length(blocks[[species]]$foldID)
# [1] 10019
# > nrow(pa_nat)
# [1] 10004
# 
# 
# esttro
# Error in `$<-.data.frame`(`*tmp*`, "foldID", value = c(2L, 3L, 1L, 5L,  )): 
#   replacement has 10086 rows, data has 10081
# 
# > length(blocks[[species]]$foldID)
# [1] 10086
# > nrow(pa_nat)
# [1] 10081
# 
# 
# plomel
# Error in `$<-.data.frame`(`*tmp*`, "foldID", value = c(1L, 1L, 4L, 3L,  )): 
#   replacement has 10109 rows, data has 10103
# 
# > length(blocks[[species]]$foldID)
# [1] 10109
# > nrow(pa_nat)
# [1] 10103
# 
# 
# psieup
# Error in `$<-.data.frame`(`*tmp*`, "foldID", value = c(5L, 5L, 4L, 1L,  )): 
#   replacement has 10277 rows, data has 10254
# 
# 
# psikra
# Error in `$<-.data.frame`(`*tmp*`, "foldID", value = c(4L, 4L, 4L, 2L,  )): 
#   replacement has 10803 rows, data has 10785
# 
# 
# species_list
# species_list[5:length(species_list)]
species_list[(grep(species, species_list)+1):length(species_list)]
