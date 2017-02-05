data <- data.exp.grouped.equal.all
path.optimisation.data <- path.optimisation.data
plot.title <- ""
filename.data_model = "data_model_all_stm.csv"
grid.ncol = 4
g.width <- 24
g.height <- 12
#### experimental data ####
data.summarise <- data %>% 
  group_by(stimulation, priming, time) %>% 
  summarise(intensity_mean = mean(intensity), 
            intensity_sd = sd(intensity))

data.summarise <- data.summarise %>% mutate(intensity_lmvn_mean = mean.lmvn(intensity_mean, intensity_sd), 
                                            intensity_lmvn_sd = sd.lmvn(intensity_mean, intensity_sd))

#### model data  ####
id <- optimisation.best[1]
data.model.best <- read.table(
  file = paste(path.optimisation.data, id, filename.data_model, sep = "/"),
  sep = ",",
  header = TRUE)  
data.model.best$id <- "optimised"

names(data.model.list) <- c("initial", "double receptors")
for( name in names(data.model.list) ){
  data.model.list[[name]]$id <- name
  data.model.best <- rbind(data.model.best, data.model.list[[name]])
} 

#### normal  ####
prm <- 0
g <- ggplot(data = data.summarise %>% filter(stimulation != 0, priming == prm), aes(x = factor(time - 5), y = intensity_mean)) + 
  geom_errorbar(aes(ymin = intensity_mean - sqrt(intensity_sd), ymax = intensity_mean + sqrt(intensity_sd))) +
  geom_point() +
  facet_grid(~stimulation) +
  geom_point(data = data.model.best %>% filter(priming == prm), aes(x = factor(time - 5), y = m.norm, color = id, group = id)) +
  geom_line(data = data.model.best %>% filter(priming == prm), aes(x = factor(time - 5), y = m.norm, color = id, group = id)) +
  ggtitle("Model comparison. Control") +
  theme_jetka() +
  xlab("time [min]") +
  ylab("intensity") +
  ylim(c(0,750))
ggsave(paste(path.optimisation, "best-model-control.pdf", sep = "/"), width = g.width, height =  g.height, useDingbats = FALSE)

prm <- 1000
g <- ggplot(data = data.summarise %>% filter(stimulation != 0, priming == prm), aes(x = factor(time - 5), y = intensity_mean)) + 
  geom_errorbar(aes(ymin = intensity_mean - sqrt(intensity_sd), ymax = intensity_mean + sqrt(intensity_sd))) +
  geom_point(shape = 23) +
  facet_grid(~stimulation) +
  geom_point(data = data.model.best %>% filter(priming == prm), aes(x = factor(time - 5), y = m.norm, color = id, group = id)) +
  geom_line(data = data.model.best %>% filter(priming == prm), aes(x = factor(time - 5), y = m.norm, color = id, group = id)) +
  ggtitle("Model comparison. Prestimulation IFNB") +
  theme_jetka() +
  xlab("time [min]") +
  ylab("intensity") +
  ylim(c(0,750))
ggsave(paste(path.optimisation, "best-model-priming.pdf", sep = "/"), width =  g.width, height =  g.height, useDingbats = FALSE)


#### logarythmic ####
prm <- 0
ggplot(data = data.summarise %>% filter(stimulation != 0, priming == prm), aes(x = factor(time - 5), y = intensity_lmvn_mean)) + 
  geom_errorbar(aes(ymin = intensity_lmvn_mean - sqrt(intensity_lmvn_sd), ymax = intensity_lmvn_mean + sqrt(intensity_lmvn_sd))) +
  geom_point() +
  facet_grid(~stimulation) +
  geom_point(data = data.model.best %>% filter(priming == prm), aes(x = factor(time - 5), y = log(m.norm), color = id, group = id)) +
  geom_line(data = data.model.best %>% filter(priming == prm), aes(x = factor(time - 5), y = log(m.norm), color = id, group = id)) +
  ggtitle("Model comparison. Control") +
  theme_jetka() +
  xlab("time [min]") +
  ylab("log(intensity)")

prm <- 1000
ggplot(data = data.summarise %>% filter(stimulation != 0, priming == prm), aes(x = factor(time - 5), y = intensity_lmvn_mean)) + 
  geom_errorbar(aes(ymin = intensity_lmvn_mean - sqrt(intensity_lmvn_sd), ymax = intensity_lmvn_mean + sqrt(intensity_lmvn_sd))) +
  geom_point() +
  facet_grid(~stimulation) +
  geom_point(data = data.model.best %>% filter(priming == prm), aes(x = factor(time - 5), y = log(m.norm), color = id, group = id)) +
  geom_line(data = data.model.best %>% filter(priming == prm), aes(x = factor(time - 5), y = log(m.norm), color = id, group = id)) +
  ggtitle("Model comparison. Prestimulation IFNB") +
  theme_jetka() + 
  xlab("time [min]") +
  ylab("log(intensity)")

#### ####
  
  pdf(file = paste(path.output, plot.title, ".pdf", sep = ""),
      width = 20,
      height = 12,
      useDingbats = FALSE)
  print(marrangeGrob(unlist(plot.list, recursive = FALSE), ncol = grid.ncol, nrow =  grid.nrow))
  dev.off()
  
  
  
  #### ####
  data <- data.exp.grouped.all
  
  
  
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
