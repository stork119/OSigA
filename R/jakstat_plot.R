### ### ###
### plot compare
### ### ###
source("R/theme_jetka.R")

compare_models <- function(
  data = data.exp.grouped,
  data.model.list = list(def = res$data.model),
  plot.title = "",
  filename,
  plot.save = TRUE,
  plot.return = TRUE
){
  data.model <- data.table(label = character(),
                           time  = numeric(),
                           m     = numeric(),
                           sd    = numeric(),
                           priming = numeric(),
                           stimulation = numeric(),
                           m.norm = numeric(),
                           sd.norm = numeric(),
                           mean.lmvn = numeric(),
                           sd.lmvn = numeric()
                           )
  for(data.model.i in 1:length(data.model.list)){
    data.model.list[[data.model.i]]$label <- labels(data.model.list)[data.model.i]
    data.model <- rbind(data.model,
                        data.model.list[[data.model.i]][, c("label" ,      "time" ,       "m"     ,      "sd",          "priming",     "stimulation", "m.norm",      "sd.norm",     "mean.lmvn" , "sd.lmvn")]
                        ### TODO ###
                        )

  }
  data.distinct <- data %>% ungroup() %>% distinct(priming, stimulation)
  data.model$intensity <- data.model$m.norm
  g.list <- list()
  if(plot.save){
    pdf(file = paste(filename, "-list", ".pdf", sep = ""),
        width = 16,
        height = 8,
        useDingbats = FALSE)
  }
  for(data.distinct.i in 1:nrow(data.distinct)){
    data.distinct.tmp <- data.distinct[data.distinct.i,]
    
    g.list[[data.distinct.i]] <- ggplot(data = data %>% filter(priming == data.distinct.tmp$priming,
                                  stimulation == data.distinct.tmp$stimulation),
           mapping = aes(x = factor(time - 5), y = intensity)) +
      geom_boxplot() +
      geom_line( data = data.model %>% filter(priming == data.distinct.tmp$priming,
                                        stimulation == data.distinct.tmp$stimulation),
                 mapping = aes(x = factor(time - 5), y = m.norm, group = factor(label), color = factor(label))) +
      geom_errorbar(
        data = data.model %>% filter(priming == data.distinct.tmp$priming,
                                     stimulation == data.distinct.tmp$stimulation),
        mapping = aes(x = factor(time - 5), ymax = m.norm + sqrt(sd.norm), ymin=m.norm - sqrt(sd.norm), color = factor(label))) +
      ggtitle(paste(plot.title, data.distinct.tmp$priming, data.distinct.tmp$stimulation)) +
      ylim(c(0,1000)) + 
      theme_jetka()
    if(plot.save){
      print(g.list[[data.distinct.i]])
    }
  }  
  if(plot.save){
    dev.off()
  }
  
  if(plot.save){
    pdf(file = paste(filename, "-grid", ".pdf", sep = ""),
        width = 20,
        height = 12,
        useDingbats = FALSE)
    print(marrangeGrob(g.list, ncol = 4, nrow = 2))
    dev.off()
  }
  
  if(plot.return){
    return(g.list[[data.distinct.i]])
  }
}

#### save results ####
save_results <- function(path.opt = paste(path.cmaes, "2017-01-26", sep ="/"),
                         par.opt = cma_es.res$xmin,
                         res.list,
                         data.exp.grouped,
                         fun_run_model = run_model_mean,
                         fun.likelihood = fun.likelihood.mvn.mean,
                         ...
){
  
  res.opt <- optimisation(fun_run_model = fun_run_model,
                          par = par.opt,
                          return.model = TRUE,
                          fun.likelihood = fun.likelihood,
                          data.exp.grouped = data.exp.grouped,
                          ...)
  
  dir.create(path.opt, recursive = TRUE)
  
  write.table(x = matrix(par.opt, ncol = 1), 
              file = paste(path.opt, "par.txt", sep ="/"),
              sep = ",",
              row.names = FALSE,
              col.names = FALSE)
  
  write.table(x = matrix(res.opt$optimisation, ncol = 1), 
              file = paste(path.opt, "optimisation.csv", sep ="/"),
              sep = ",",
              row.names = FALSE,
              col.names = FALSE)
  write.table(x = res.opt$data.model,
              file = paste(path.opt, "data_model.csv", sep ="/"),
              sep = ",",
              row.names = FALSE,
              col.names = TRUE)
  
  g <- compare_models(
    data = data.exp.grouped,
    data.model.list = c(list(opt = res.opt$data.model),res.list),
    plot.title = "Model compares",
    filename = paste(path.opt, "compare", sep ="/"),
    plot.save = TRUE,
    plot.return = TRUE
  )
  return(res.opt)
}

