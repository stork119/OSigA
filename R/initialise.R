###
### initialise 
###

#mesh.exp <- unique(data.exp$time)
#data.exp.grouped <- data.exp.grouped %>% filter(stimulation == 0.1 | stimulation == 1)

path.parameters <- "resources/input/"
par.def <- scan(file = paste(path.parameters, "par.txt", sep = ""))
varscale <- 0.15

tmesh <- seq(from = 0, to = 100, by = 1)
tmesh.list <- which(tmesh %in% (unique(data.list$data.exp$time)))

stimulation.list <- (data.list$data.exp %>% 
                       dplyr::distinct(stimulation) %>% 
                       dplyr::filter(stimulation != 0))$stimulation

background     <- mean(data.list$data.exp[data.list$data.exp$time == 0,]$intensity)
background.var <- var(data.list$data.exp[data.list$data.exp$time == 0,]$intensity)


variables <- rep(0.0, times = 629)
variables.priming <- rep(0.0, times = 629)
variables[1:17] <- scan(file = paste(path.parameters, "var.txt", sep = ""))
variables.priming[1:17] <- scan(file = paste(path.parameters, "var-priming.txt", sep = ""))
variables[44:52] <- scan(file = paste(path.parameters, "var-extrinsic.txt", sep = ""))
variables.priming[44:52] <- variables[44:52]
variables[27:43] <- varscale*(variables[1:17]^2)
variables.priming[27:43] <- varscale*(variables.priming[1:17]^2)


#### ####
# data.model.list <- list()
# path.single <- paste(path.output, "single", sep = "/")
# data.model.list$single <- read.table(file = paste(path.single, "data_model.csv", sep = "/"),
#                                      sep = ",", header = TRUE)
# data.model.list$single$time <-  data.model.list$single$time - 5
# path.receptors <- paste(path.output, "receptors", sep = "/")
# data.model.list$double <- read.table(file = paste(path.receptors, "data_model.csv", sep = "/"),
#                                      sep = ",", header = TRUE)
# data.model.list$double$time <-  data.model.list$double$time - 5
#### ####