# #### ####
# 
# ids <- list.dirs(path.optimisation, full.names = FALSE)
# ids <- ids[which(!is.na(as.numeric(ids)))]
# 
# registerDoParallel(no_cores)
# foreach( id = ids[which(!ids %in% optimisation.table$id)] ) %dopar%{
#   try({
#     
#     data.model <-  read.table(
#       file = paste(path.optimisation, id, "data_model.csv", sep = "/"),
#       sep = ",",
#       header = TRUE)
# 
#     result <- sapply(fun.likelihood.list,
#                      function(fun.likelihood){
#                        sum( likelihood(
#                          fun.likelihood = fun.likelihood,
#                          data.model = data.model,
#                          data.exp.grouped = data.exp.grouped))
#                      })
# 
#     write.table(x = matrix(result, nrow = 1), paste(path.optimisation, id, "optimisation.csv", sep = "/"), sep = ",", row.names = FALSE, col.names = FALSE)
# 
#   })
# }
# stopImplicitCluster()
# 

#### ####
#ids <- list.dirs(path.optimisation.data, full.names = FALSE)
source("R/sum_data_initialise.R")
ids.opt <- ids[-which(ids %in% par.list.ids)]
ids.opt  <- ids.opt [which(!is.na(as.numeric(ids.opt )))]


for( id  in ids.opt [which(!ids.opt  %in% optimisation.table$id)] ){
  try({
        optimisation.table <- rbind(optimisation.table,
                                read_optimisation(path = paste(path.optimisation.data, id, sep = "/"),
                                                  id = id, 
                                                  names = names(fun.likelihood.list)))
  })
}

#### ####
fun.likelihood.name <- "sd_data"
optimisation.table <- optimisation.table[order(as.numeric(optimisation.table[,fun.likelihood.name]) ),]
optimisation.best <- c(optimisation.table[order(as.numeric(optimisation.table[,fun.likelihood.name])[1:6]),]$id, "single", "receptors")
no_cores <- 4
#stimulation.list <- (data.exp.grouped %>% ungroup() %>% distinct(stimulation))$stimulation
# data.exp.grouped.all <- data.list$data.exp.norm %>% group_by(priming, stimulation, time) %>% mutate(intensity_sd = var(intensity))
# data.exp.grouped.equal.all <- get_equal_data(data.exp.grouped.all)
# registerDoParallel(no_cores)
# optimisation.table.all <- foreach( i = 1:length(optimisation.best), .combine = rbind) %dopar% {
#   id <- optimisation.best[i]
#   try({
#     parameters <- scan(paste(path.optimisation.data, id, "par.txt", sep = "/"))
#     variables <- scan(paste(path.optimisation.data, id, "var.txt", sep = "/"))
#     variables.priming <- scan(paste(path.optimisation.data, id, "var-priming.txt", sep = "/"))
#     filename <- paste(path.optimisation.data, id, "data_model_all_stm.csv", sep = "/")
#     model.simulation <- list()
#     if(file.exists(filename)){
#       model.simulation$data.model <- read.table(
#                   file = filename,
#                   sep = ",",
#                   header = TRUE)
#     } else {
#       model.simulation <- do.call(run_model_mean,
#                                 list(parameters = parameters,
#                                      variables = variables,
#                                      variables.priming = variables.priming,
#                                      tmesh = tmesh,
#                                      tmesh.list = tmesh.list,
#                                      stimulation.list = stimulation.list.all,
#                                      background = background))
#       write.table(x = model.simulation$data.model, 
#                 file = filename,
#                 sep = ",",
#                 row.names = FALSE,
#                 col.names = TRUE)
#     }
#     result <- sapply(fun.likelihood.list,
#                      function(fun.likelihood){
#                        sum( likelihood(
#                          fun.likelihood = fun.likelihood,
#                          data.model = model.simulation$data.model,
#                          data.exp.grouped = data.exp.grouped.optimisation,
#                          data.exp.summarise  = data.exp.summarise.optimisation))
#                      })
#     
#     write.table(x = matrix(result, nrow = 1), 
#                 paste(path.optimisation.data, id, "optimisation_all.csv", sep = "/"),
#                 sep = ",", row.names = FALSE, col.names = FALSE)
#     return( matrix(c(id,result), nrow = 1))
#   })
# }
# stopImplicitCluster()

registerDoParallel(no_cores)
data.model.list <- foreach( i = 1:length(optimisation.table$sd_data)) %dopar% {
  id <- optimisation.table$id[i]
  try({
    filename <- paste(path.optimisation.data, id, "data_model.csv", sep = "/")
    if(file.exists(filename)){
      data.model <- read.table(
                  file = filename,
                  sep = ",",
                  header = TRUE)
      data.model$type <- id
    return( data.model )
    }
  })
}
stopImplicitCluster()

data.model <- do.call(rbind, data.model.list)

optimisation.table.results <- data.model %>% dplyr::group_by(type) %>% summarise(likelihood = sum(likelihood)) %>% mutate(id = type) %>% arrange(likelihood)
#colnames(optimisation.table.all) <- colnames(optimisation.table)
write.table(file = paste(path.optimisation.results, "optimisation_ranking_all.csv", sep = ""),
            x = optimisation.table.results,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)



registerDoParallel(no_cores)
data.parameters.list <- foreach( i = 1:length(optimisation.table$sd_data)) %dopar% {
  id <- optimisation.table$id[i]
  try({
    filename <- paste(path.optimisation.data, id, "par_exp.txt", sep = "/")
    if(file.exists(filename)){
      parameters.id <- scan(file = filename)
      data.parameters <- matrix(parameters.id, nrow = 1)
      colnames(data.parameters) <- paste("p", 1:ncol(data.parameters), sep = "")
      return( data.parameters %>% data.frame() %>% mutate(id = id) )
    }
  })
}
stopImplicitCluster()
data.parameters <- do.call(rbind, data.parameters.list)

data.parameters.opt <- data.parameters[,par.optimised]
data.parameters.opt$id <- data.parameters$id

data.parameters.opt <- data.parameters.opt %>% left_join(optimisation.table.results)
ggplot(data.parameters.opt %>% 
         melt(id.vares = id) %>% 
         left_join(optimisation.table.results) %>% 
         filter(!(id %in% c("single", "receptors"))),
       aes(y = value, x = variable, group = id, colour = likelihood ) )+ geom_point() + geom_line()

parameters.table.all <- foreach( i = 1:nrow(optimisation.table), .combine = rbind) %do% {
  id <- optimisation.table$id[i]
  try({
    parameters <- scan(paste(path.optimisation.data, id, "par.txt", sep = "/"))
    return( matrix(c(id,parameters), nrow = 1))
  })
}
colnames(parameters.table.all) <- c("id", 1:10)
write.table(file = paste(path.optimisation, "parameters_ranking.csv", sep = ""),
            x = parameters.table.all,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
write.table(file = paste(path.optimisation, "optimisation_ranking.csv", sep = ""),
            x = optimisation.table,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)

plot_results(data = data.exp.grouped.equal.all,
             path.analysis = path.optimisation.data,
             path.output = path.optimisation,
             data.model.best.list = data.model.list,
             optimisation.best = optimisation.best[1:5],#st.likelihood[order(st.likelihood$priming),]$id[1:10],
             plot.title = "all_fullmodel", 
             filename.data_model = "data_model_all_stm.csv",
             grid.ncol = 4)


plot_results(data = data.exp.grouped.equal.all,
             path.analysis = path.optimisation.data,
             path.output = path.optimisation,
             data.model.best.list = data.model.list,
             optimisation.best = optimisation.best[1:5],#st.likelihood[order(st.likelihood$priming),]$id[1:10],
             plot.title = "all", 
             filename.data_model = "data_model.csv",
             grid.ncol = 4)


plot_results(data = data.exp.grouped.equal.all,
             path.analysis = path.optimisation.data,
             path.output = path.optimisation,
             data.model.best.list = data.model.list,
             optimisation.best = st.likelihood[order(st.likelihood$priming),]$id[1:10],
             plot.title = "priming", 
             filename.data_model = "data_model.csv",
             grid.ncol = 4)




#### ####
# plot_results(data = data.exp.grouped,
#              path.optimisation = path.optimisation,
#              data.model.best.list = data.model.list,
#              optimisation.best =  optimisation.table[order(optimisation.table$mvn.mean)[1:5],]$id,
#              plot.title = "mvn_mean_all", 
#              filename.data_model = "data_model_all_stm.csv",
#              grid.ncol = 4)
# 
# plot_results(data = data.exp.grouped,
#              path.optimisation = path.optimisation,
#              data.model.best.list = data.model.list,
#              optimisation.best = optimisation.table[order(optimisation.table$lmvn)[1:5],]$id,
#              plot.title = "lmvn")
# 
# 
# parameters.matrix <- matrix(par.def, ncol = 10, nrow = 1)
# optimisation.best <- optimisation.table[order(optimisation.table$mvn.mean)[1:5],]$id
# for( i in 1:length(optimisation.best)){
#   id <- optimisation.best[i]
#   try({
#     parameters.matrix <- rbind(parameters.matrix, scan(paste(path.optimisation, id, "par.txt", sep = "/")))
# 
#   })
# }
