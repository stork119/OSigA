#### ####

ids <- list.dirs(path.optimisation, full.names = FALSE)
ids <- ids[which(!is.na(as.numeric(ids)))]

for( id in ids[which(!ids %in% optimisation.table$id)] ) {
  try({
  optimisation.table <- rbind(optimisation.table,
                              read_optimisation(path = paste(path.optimisation, id, sep = "/"),
                                                id = id))
  })
}


optimisation.best <- optimisation.table[order(optimisation.table$mvn.mean)[1:5],]$id
#stimulation.list <- (data.exp.grouped %>% ungroup() %>% distinct(stimulation))$stimulation
for( i in 1:length(optimisation.best)){
  id <- optimisation.best[i]
  try({
    parameters <- scan(paste(path.optimisation, id, "par.txt", sep = "/"))
    model.simulation <- do.call(run_model,
                                list(parameters = parameters,
                                     variables = variables,
                                     variables.priming = variables.priming,
                                     tmesh = tmesh,
                                     tmesh.list = tmesh.list,
                                     stimulation.list = stimulation.list,
                                     background = background))
    write.table(x = model.simulation$data.model, 
                file = paste(path.optimisation, id, "data_model_all_stm.csv", sep = "/"),
                sep = ",",
                row.names = FALSE,
                col.names = TRUE)
    
  })
}

plot_results(data = data.exp.grouped.all,
             path.optimisation = path.optimisation,
             data.model.best.list = data.model.list,
             optimisation.best =  optimisation.best,
             plot.title = "mvn.mean_all", 
             filename.data_model = "data_model_all_stm.csv",
             grid.ncol = 4)


optimisation.table <- optimisation.table[order(optimisation.table$mvn),]
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
