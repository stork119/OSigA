### ###
### 
### ###

id <- "204"
analyse_optimisation <- function(path.optimisation.data,
                                 id,
                                 return.rds = FALSE,
                                 path.optimisation.results){
  
  path.optimisation.results.data <- paste(path.optimisation.results, "data", id, sep = "/")
  dir.create(path = path.optimisation.results.data, recursive = TRUE, showWarnings = FALSE)
  optimisation.analyse <- list()
  optimisation.plots <- list()
  optimisation.list <- readRDS(file = paste(path.optimisation.data, id, "optimisation.rds", sep = "/"))
  if(return.rds){
    optimisation.analyse[["rds"]]  <- optimisation.list
  }
  optimisation.plots[["sigma"]] <-
    ggplot(data = data.frame(step  = 1:length(optimisation.list$diagnostic$sigma),
                             sigma = optimisation.list$diagnostic$sigma),
           mapping = aes(x = step, y = sigma)) + 
    geom_point() + 
    do.call(theme_jetka, plot.args) + 
    ggtitle("sigma")


  optimisation.plots[["likelihood"]] <-
    ggplot(data = data.frame(step  = 1:length(optimisation.list$diagnostic$sigma),
                             values = sapply(1:nrow(optimisation.list$diagnostic$value),
                                             function(step){min(optimisation.list$diagnostic$value[step,])})),
           mapping = aes(x = step, y = log(values))) + 
    geom_point()+ 
    do.call(theme_jetka, plot.args) + 
    ggtitle("Likelihood")

  data.optimisation.list  <- foreach(step = 1:dim(optimisation.list$diagnostic$pop)[3]) %do% {
    return(data.frame(step = step,
                      parameter = paste("p", par.optimised, sep = ""),
                      value = optimisation.list$diagnostic$pop[,which.min(optimisation.list$diagnostic$value[step,]),step],
                      likelihood = min(optimisation.list$diagnostic$value[step,]),
                      
                      sigma = optimisation.list$diagnostic$sigma[step]
                      ))
  }
  data.optimisation <- do.call(what = rbind, args = data.optimisation.list) %>% data.table()

  optimisation.plots[["steps"]] <- 
    ggplot( 
      data.optimisation,
      aes(y = value, x = parameter, group = step, colour = step ) ) +
    geom_point() +
    geom_line() + 
    do.call(theme_jetka, plot.args) +
    ggtitle("Values of parameters in steps")
  
  
  optimisation.plots[["steps_likelihood"]] <- ggplot( 
    data.optimisation,
    aes(y = value, x = -log(likelihood), group = parameter, colour = parameter ) ) +
    geom_point() +
    geom_line() + 
    do.call(theme_jetka, plot.args) +
    ggtitle("Likelihood vs parameters values")
  
  
  do.call(what = ggsave, 
          args = append(plot.args.ggsave, 
                        list(filename = paste(path.optimisation.results.data, "optimisation.pdf", sep = "/"),
                             plot = marrangeGrob(grobs = optimisation.plots, ncol = 1, nrow = 1))))
  optimisation.analyse[["plots"]] <- optimisation.plots
  return(optimisation.analyse)
}