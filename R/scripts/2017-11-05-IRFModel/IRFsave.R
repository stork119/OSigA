### ###
### saving results 
### ###


saveResults <- function(
  model.type = "irf", #"pstat"
  irfmodel.path.list,
  output.path = NULL,
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
  if(is.null(output.path)){
    output.path <- irfmodel.path.list$output.path
  }
  dir.create(output.path, 
             recursive = TRUE)
  
  if(!is.null(g.list)){
    do.call(what = ggsave,
            args = append(plot.args.ggsave,
                          list(filename = paste(output.path, 
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
    file = paste(output.path, 
                 paste("IRFmodel-", 
                       model.type,
                       ".RDS", 
                       sep = ""),
                 sep = "/"),
    object = results)   
  return()
}