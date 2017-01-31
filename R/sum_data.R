#### ####

ids <- list.dirs(path.optimisation, full.names = FALSE)
ids <- ids[which(!is.na(as.numeric(ids)))]

registerDoParallel(no_cores)
foreach( id = ids[which(!ids %in% optimisation.table$id)] ) %dopar%{
  try({
    
    data.model <-  read.table(
      file = paste(path.optimisation, id, "data_model.csv", sep = "/"),
      sep = ",",
      header = TRUE)

    result <- sapply(fun.likelihood.list,
                     function(fun.likelihood){
                       sum( likelihood(
                         fun.likelihood = fun.likelihood,
                         data.model = data.model,
                         data.exp.grouped = data.exp.grouped))
                     })

    write.table(x = matrix(result, nrow = 1), paste(path.optimisation, id, "optimisation.csv", sep = "/"), sep = ",", row.names = FALSE, col.names = FALSE)

  })
}
stopImplicitCluster()


#### ####
ids <- list.dirs(path.optimisation, full.names = FALSE)
ids <- ids[which(!is.na(as.numeric(ids)))]

for( id in ids[which(!ids %in% optimisation.table$id)] ){
  try({
        optimisation.table <- rbind(optimisation.table,
                                read_optimisation(path = paste(path.optimisation, id, sep = "/"),
                                                  id = id))
  })
}

#### ####
optimisation.best <- c(optimisation.table[order(optimisation.table$mvn.sd_const)[1:5],]$id, "single", "receptors")
#stimulation.list <- (data.exp.grouped %>% ungroup() %>% distinct(stimulation))$stimulation
data.exp.grouped.all <- data.exp.grouped.all %>% group_by(priming, stimulation, time) %>% mutate(intensity_sd = var(intensity))
data.exp.grouped.equal.all <- get_equal_data(data.exp.grouped.all)
registerDoParallel(no_cores)
optimisation.table.all <- foreach( i = 1:length(optimisation.best), .combine = rbind) %dopar% {
  id <- optimisation.best[i]
  try({
    parameters <- scan(paste(path.optimisation, id, "par.txt", sep = "/"))
    variables <- scan(paste(path.optimisation, id, "var.txt", sep = "/"))
    variables.priming <- scan(paste(path.optimisation, id, "var-priming.txt", sep = "/"))
    filename <- paste(path.optimisation, id, "data_model_all_stm.csv", sep = "/")
    model.simulation <- list()
    if(file.exists(filename)){
      model.simulation$data.model <- read.table(
                  file = filename,
                  sep = ",",
                  header = TRUE) 
    } else {
      model.simulation <- do.call(run_model,
                                list(parameters = parameters,
                                     variables = variables,
                                     variables.priming = variables.priming,
                                     tmesh = tmesh,
                                     tmesh.list = tmesh.list,
                                     stimulation.list = stimulation.list.all,
                                     background = background))
      write.table(x = model.simulation$data.model, 
                file = filename,
                sep = ",",
                row.names = FALSE,
                col.names = TRUE)
    }
    result <- sapply(fun.likelihood.list,
                     function(fun.likelihood){
                       sum( likelihood(
                         fun.likelihood = fun.likelihood,
                         data.model = model.simulation$data.model,
                         data.exp.grouped = data.exp.grouped.equal.all))
                     })
    
    write.table(x = matrix(result, nrow = 1), 
                paste(path.optimisation, id, "optimisation_all.csv", sep = "/"),
                sep = ",", row.names = FALSE, col.names = FALSE)
    return( matrix(c(id,result), nrow = 1))
  })
}
stopImplicitCluster()
colnames(optimisation.table.all) <- colnames(optimisation.table)
write.table(file = paste(path.optimisation, "optimisation_ranking_all.csv", sep = ""),
            x = optimisation.table.all,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)

parameters.table.all <- foreach( i = 1:length(optimisation.best), .combine = rbind) %do% {
  id <- optimisation.best[i]
  try({
    parameters <- scan(paste(path.optimisation, id, "par.txt", sep = "/"))
    return( matrix(c(id,parameters), nrow = 1))
  })
}
colnames(parameters.table.all) <- c("id", 1:10)
write.table(file = paste(path.optimisation, "parameters_ranking__mvn_mean_all.csv", sep = ""),
            x = parameters.table.all,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)


plot_results(data = data.exp.grouped.equal.all,
             path.analysis = path.optimisation,
             data.model.best.list = data.model.list,
             optimisation.best =  optimisation.best[1:5],
             plot.title = "mvn.sd_const_all", 
             filename.data_model = "data_model_all_stm.csv",
             grid.ncol = 4)


optimisation.table <- optimisation.table[order(optimisation.table$mvn.sd_const),]
write.table(file = paste(path.optimisation, "optimisation_ranking.csv", sep = ""),
            x = optimisation.table,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)

plot_results(data = data.exp.grouped,
             path.optimisation = path.optimisation,
             data.model.best.list = data.model.list,
             optimisation.best =  optimisation.table[order(optimisation.table$mvn.mean)[1:5],]$id,
             plot.title = "mvn_mean_all", 
             filename.data_model = "data_model_all_stm.csv",
             grid.ncol = 4)

plot_results(data = data.exp.grouped,
             path.optimisation = path.optimisation,
             data.model.best.list = data.model.list,
             optimisation.best = optimisation.table[order(optimisation.table$lmvn)[1:5],]$id,
             plot.title = "lmvn")

#### ####
parameters.matrix <- matrix(par.def, ncol = 10, nrow = 1)
optimisation.best <- optimisation.table[order(optimisation.table$mvn.mean)[1:5],]$id
for( i in 1:length(optimisation.best)){
  id <- optimisation.best[i]
  try({
    parameters.matrix <- rbind(parameters.matrix, scan(paste(path.optimisation, id, "par.txt", sep = "/")))

  })
}
