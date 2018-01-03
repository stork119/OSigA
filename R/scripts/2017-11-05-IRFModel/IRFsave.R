### ###
### saving results 
### ###


saveResults <- function(
  model.type = "irf", #"pstat"
  irfmodel.path.list,
  optimisation.res,
  likelihood, 
  g.list = NULL,
  results = list(),
  par, 
  ranges,
  stimulations,
  stopfitness,
  fun.optimisation,
  maxit,
  ...
){
  irfmodel.path.list$output.path <-
    paste(irfmodel.path.list$output.dir,
          irfmodel.path.list$optimisation.id, sep = "/")
  
  dir.create(irfmodel.path.list$output.path, 
             recursive = TRUE)
  
  if(!is.null(g.list)){
    do.call(what = ggsave,
            args = append(plot.args.ggsave,
                          list(filename = paste(irfmodel.path.list$output.path, 
                                                paste("IRFmodel-", 
                                                      model.type,
                                                      ".pdf", 
                                                      sep = ""),
                                                sep = "/"),
                               plot = marrangeGrob(grobs = g.list, ncol = 1, nrow = 1))))
  }
  
  results <- append(results, 
                    list(
                      list(
                        optimisation = optimisation.res,
                        likelihood = likelihood,
                        par = par, 
                        ranges = ranges,
                        stimulations = stimulations,
                        stopfitness = stopfitness,
                        fun.optimisation = fun.optimisation,
                        maxit = maxit,
                        plot = g.list
                      )
                    )
  )
  
  saveRDS(
    file = paste(irfmodel.path.list$output.path, 
                 paste("IRFmodel-", 
                       model.type,
                       ".RDS", 
                       sep = ""),
                 sep = "/"),
    object = results)   
  return(irfmodel.path.list)
}