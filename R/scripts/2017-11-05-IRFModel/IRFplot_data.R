### ###
### IRFplotdata
### ###

glist.data <- list()

plot.title <- "pSTAT data"

glist.data[["pSTAT"]] <- 
  ggplot(irfmodel.data.list$pSTATsum,
       mapping = 
         aes(x = factor(stimulation), 
             y = logresponse, 
             ymin = logresponse - logresponse.sd,
             ymax = logresponse + logresponse.sd,
             group = factor(stimulation))) +
  geom_errorbar() + 
  geom_point() +
  do.call(theme_jetka, args = plot.args) +
  ggtitle(paste(plot.title))


plot.title <- "IRF data"

glist.data[["IRF"]] <- 
  ggplot(irfmodel.data.list$irfsum,
         mapping = 
           aes(x = factor(stimulation), 
               y = logresponse, 
               ymin = logresponse - logresponse.sd,
               ymax = logresponse + logresponse.sd,
               group = factor(stimulation))) +
  geom_errorbar() + 
  geom_point() +
  do.call(theme_jetka, args = plot.args) +
  ggtitle(paste(plot.title))