### ### ###
### plot compare
### ### ###
library(ggplot2)

compare_models <- function(
  data = data.exp.grouped,
  data.model.list = list(def = res$data.model)
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
  data.distinct <- data %>% distinct(priming, stimulation)
  data.model$intensity <- data.model$m.norm
  for(data.distinct.i in 1:nrow(data.distinct)){
    data.distinct.tmp <- data.distinct[data.distinct.i,]
    
    g <- ggplot(data = data %>% filter(priming == data.distinct.tmp$priming,
                                  stimulation == data.distinct.tmp$stimulation),
           mapping = aes(x = factor(time), y = intensity)) +
      geom_boxplot() +
      geom_line( data = data.model %>% filter(priming == data.distinct.tmp$priming,
                                        stimulation == data.distinct.tmp$stimulation),
                 mapping = aes(x = factor(time), y = m.norm, group = factor(label), color = factor(label))) +
      geom_errorbar(
        data = data.model %>% filter(priming == data.distinct.tmp$priming,
                                     stimulation == data.distinct.tmp$stimulation),
        mapping = aes(x = factor(time), ymax = m.norm + sqrt(sd.norm), ymin=m.norm - sqrt(sd.norm)))
    
  }    
  
}
