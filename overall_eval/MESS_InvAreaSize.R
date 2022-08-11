# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
setwd(paste0("../../Modeling_bird_invasions/PA_modelling/inv_region_tables/"))
getwd()

# LOAD PACKAGES ####


###################################################################################################################################
###READ IN DATA####################################################################################################################
###################################################################################################################################
temp <- list.files(pattern="*.csv")
mess.ir <- lapply(temp, read.csv)

###################################################################################################################################
###Analyses########################################################################################################################
###################################################################################################################################  
summary.data <- data.frame(matrix(vector(), 0, 2))
colnames(summary.data) <- c("eu.mess.ratio","ir.mess.ratio")

for (i in 1:20){
mess.ir[[i]]
  eu <- mess.ir[[i]]
  eu.mess <- mess.ir[[i]][mess.ir[[i]]$analog ==1,]
  
  ir <- mess.ir[[i]][mess.ir[[i]]$inv_range ==1,]
  ir.mess <- ir[ir$analog==1,]
  
  eu.mess.ratio <- nrow(eu.mess)/nrow(eu)
  ir.mess.ratio <- nrow(ir.mess)/nrow(ir)
  temp.ratios <- cbind(eu.mess.ratio,ir.mess.ratio)
  summary.data <- rbind(summary.data,temp.ratios)
}
 
mean(summary.data$eu.mess.ratio)
sd(summary.data$eu.mess.ratio)

mean(summary.data$ir.mess.ratio)
sd(summary.data$ir.mess.ratio)