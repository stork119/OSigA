### ### ###
### plot compare
### ### ###
source("R/graphics/libraries.R")

compare_models <- function(
  data,
  data.model.list = list(def = res$data.model),
  plot.title = "",
  filename,
  plot.save = TRUE,
  plot.return = TRUE,
  grid.nrow = 2,
  grid.ncol = max(as.numeric(unique(data$stimulation)/grid.nrow),1),
  ...
){
  
  data.summarise <- data %>% 
    group_by(stimulation, priming, time) %>% 
    summarise(intensity_mean = mean(intensity), 
              intensity_sd = sd(intensity))
  
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
      geom_point(data = data.summarise %>% 
                   filter(priming == data.distinct.tmp$priming,
                          stimulation == data.distinct.tmp$stimulation),
                 aes(x = factor(time - 5), y = intensity_mean), shape = 23) +
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
    print(marrangeGrob(g.list, ncol = grid.ncol, nrow =  grid.nrow))
    dev.off()
  }
  
  if(plot.return){
    return(g.list)
  }
}

#### save results ####
save_results <- function(path.opt,
                         data.model.opt,
                         optimisation.opt,
                         optimisation.opt.colnames = FALSE,
                         par.opt,
                         par.exp.opt = par.opt,
                         res.list,
                         data.exp.grouped,
                         error = FALSE,
                         variables = matrix(),
                         variables.priming = matrix(),
                         ...
){
  
  dir.create(path.opt, recursive = TRUE)
  
  if(error){
    write.table(x = matrix(1, ncol = 1), 
              file = paste(path.opt, "error.txt", sep ="/"),
              sep = ",",
              row.names = FALSE,
              col.names = FALSE)
  }

  write.table(x = matrix(variables, ncol = 1), 
              file = paste(path.opt, "var.txt", sep ="/"),
              sep = ",",
              row.names = FALSE,
              col.names = FALSE)
  
  write.table(x = matrix(variables.priming, ncol = 1), 
              file = paste(path.opt, "var-priming.txt", sep ="/"),
              sep = ",",
              row.names = FALSE,
              col.names = FALSE)
    
  write.table(x = matrix(par.opt, ncol = 1), 
              file = paste(path.opt, "par.txt", sep ="/"),
              sep = ",",
              row.names = FALSE,
              col.names = FALSE)
  
  write.table(x = matrix(par.exp.opt, ncol = 1), 
              file = paste(path.opt, "par_exp.txt", sep ="/"),
              sep = ",",
              row.names = FALSE,
              col.names = FALSE)
  
  write.table(x = matrix(optimisation.opt, nrow = 1), 
              file = paste(path.opt, "optimisation.csv", sep ="/"),
              sep = ",",
              row.names = FALSE,
              col.names = optimisation.opt.colnames)
  
  write.table(x = data.model.opt,
              file = paste(path.opt, "data_model.csv", sep ="/"),
              sep = ",",
              row.names = FALSE,
              col.names = TRUE)
  
  g <- compare_models(
    data = data.exp.grouped,
    data.model.list = c(list(optimisation = data.model.opt),
                        res.list),
    plot.title = "Model compares",
    filename = paste(path.opt, "compare", sep ="/"),
    plot.save = TRUE,
    plot.return = TRUE,
    ...
  )
}

