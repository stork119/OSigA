### ###
### IRF channel capacity
### ###


library(e1071)
library(CapacityLogReg)
source("R/scripts/2017-11-05-IRFModel/IRFlibrary.R")

#### read model ####
#irfmodel.path.list$optimisation.id <- "2017-12-14-summarise-id"
irfmodel.path.list$output.path <-
  paste(irfmodel.path.list$output.dir,
        irfmodel.path.list$optimisation.id, sep = "/")
results <- readRDS(
  file = paste(irfmodel.path.list$output.path, "IRFmodel.RDS", sep = "/"))

par <- results$par#[ranges.opt]
par.new <- par
ranges.factor <- results$ranges.factor
ranges.base <- results$ranges.base
ranges.opt <- results$ranges.opt
params <- ranges.factor
params[ranges.opt] <- ranges.factor[ranges.opt]*ranges.base[ranges.opt]^par.new

stimulations <- results$stimulations

data.model <- model_fun_stm_params(
  stimulations = stimulations,
  params = params)

data <- data.model %>% dplyr::left_join(data.raw.sum, by = "stimulation", suffix = c(".model", ".data"))
#data <- data[data$stimulation != 0.01, ]


# data <-
#   data %>%
#   dplyr::mutate(pstat.mean = pstat.model.model,
#                 pstat.sd   =   pstat.sd.data,
#                 irf.mean   =  irf.model.model,
#                 irf.sd     = irf.sd.data) %>%
#   dplyr::select(stimulation,
#                 pstat.mean,
#                 pstat.sd,
#                 irf.mean,
#                 irf.sd)

data <-
  data %>%
  dplyr::mutate(pstat.mean = pstat.model,
                pstat.sd   =   pstat.sd.data,
                irf.mean   =  irf.model,
                irf.sd     = irf.sd.data) %>%
  dplyr::select(stimulation,
                pstat.mean,
                pstat.sd,
                irf.mean,
                irf.sd)

#### model sampling ####
model.sampling <- function(
  data,
  n.sample
){
  sample.list <- list()
  sample.list <- foreach(i = 1:nrow(data)) %do% {
    data.frame(signal   = data[i,]$stimulation,
               response.pstat = 
                 exp(
                   rnorm(n = n.sample,
                         mean = data[i,]$pstat.mean,
                         sd =  data[i,]$pstat.mean*data[i,]$pstat.sd)
                 ),
               response.irf = 
                 exp(
                   rnorm(n = n.sample,
                         mean = data[i,]$irf.mean,
                         sd   = data[i,]$pstat.mean*data[i,]$irf.sd)
                 ))
  }
  sample.df <- do.call(rbind, sample.list)
  return(sample.df)
}
#### ####
no_cores <- 6
n.sample <- 1000
rep.sample <- 25

# cc.id <- paste("_model", "n", n.sample, "rep", rep.sample, sep = "_")
# sd.list.pstat <- seq(from = 0, to = 0.05, by = 0.001)
# sd.list.irf <- seq(from = 0, to = 0.05, by = 0.001)

###sd.list mowi jaki procent sredniej jest odchylenie standardowe
cc.id <- paste("_data", "n", n.sample, "rep", rep.sample, sep = "_")
sd.list.pstat <- seq(from = 0, to = 0.5, by = 0.01)#*(params[c(9)]^2)
sd.list.irf <- seq(from = 0, to = 0.5, by = 0.01)#*(params[c(10)]^2)

sd.list <- expand.grid(irf = sd.list.irf, pstat = sd.list.pstat)

# sd.list <- data.frame(irf   = sd.list.irf,
#                       pstat = sd.list.pstat)

registerDoParallel(no_cores)
cc.list <- foreach( sd.i = 1:nrow(sd.list)) %dopar% {
  
  sd.irf <- sd.list[sd.i,]$irf
  data$irf.sd <- sd.irf
  sd.pstat <- sd.list[sd.i,]$pstat
  data$pstat.sd <- sd.pstat
  
  cc.sample.list <- foreach( sample.i = 1:rep.sample) %do% {
    
    sample.df <- model.sampling(
      data = data,
      n.sample = n.sample)
    
    col_signal   <- "signal"
    col_response <- "response.pstat"
    
    cc.df <- data.frame(
      sd.pstat =  sd.pstat,
      sd.irf   =  sd.irf,
      id = sample.i)
    
    cc.output.pstat <- 
      capacity_logreg_main( 
        sample.df,
        graphs = FALSE,
        signal = col_signal,
        response = col_response,
        model_out = FALSE,
        output_path = paste(irfmodel.path.list$output.path, "/", sep = "")
      )
    cc.df$pstat <- cc.output.pstat$cc
    
    col_response <- "response.irf"
    cc.output.irf <- 
      capacity_logreg_main(
        sample.df,
        graphs = FALSE,
        signal = col_signal,
        response = col_response,
        model_out = FALSE,
        output_path = paste(irfmodel.path.list$output.path, "/", sep = "")
      )
    cc.df$irf <- cc.output.irf$cc
    return(cc.df)
  }
  return(do.call(rbind, cc.sample.list))
}
stopImplicitCluster()

cc.df <- do.call(rbind, cc.list)

#### plotting chhannel capacity ####
g.cc.list <- list()

### together
#cc.df$num <- 1:nrow(cc.df)
cc.df.melt <- reshape2::melt(cc.df, 
                             id.vars = c("id", "sd.irf", "sd.pstat"),
                             measure.vars = c("irf", "pstat"))
cc.df.melt.sum <- cc.df.melt %>% 
  dplyr::group_by(variable, sd.irf, sd.pstat) %>%
  dplyr::summarise(mean = mean(value),
                   sd = sd(value))

g.cc.list[["all_together"]] <- 
  ggplot(cc.df.melt.sum, 
         aes_string( y = "mean",
                     ymin = "mean - sd",
                     ymax = "mean + sd",
                     color = "variable",
                     x = "sd.irf")) +
  ylim(c(0,3)) +
  ylab("Channel capacity") + 
  do.call(theme_jetka, args = plot.args)


if(rep.sample > 5){
  g.cc.list[["all_together"]]  <- g.cc.list[["all_together"]] + geom_errorbar()
} else {
  g.cc.list[["all_together"]]  <- g.cc.list[["all_together"]] + geom_point()
}

### pstat
col_response <- "pstat"
col_noise <- "sd.pstat"
cc.df.sum <- cc.df %>% 
  dplyr::group_by_(col_noise) %>%
  dplyr::summarise_("mean" = paste("mean(", col_response, ")"),
                   "sd" = paste("sd(", col_response, ")"))

g.cc.list[[col_response]] <- 
  ggplot(cc.df.sum, 
       aes_string( y = "mean",
                   ymin = "mean - sd",
                   ymax = "mean + sd",
                   x = col_noise)) +
  ggtitle(col_response) +
  ylim(c(0,3)) +
  ylab("Channel capacity") +
  do.call(theme_jetka, args = plot.args)

if(rep.sample > 5){
  g.cc.list[[col_response]]  <- g.cc.list[[col_response]] + geom_errorbar()
} else {
  g.cc.list[[col_response]]  <- g.cc.list[[col_response]] + geom_point()
}

### irf
col_response <- "irf"
col_noise <- "sd.irf"
cc.df.sum <- cc.df %>% 
  dplyr::group_by_(col_noise) %>%
  dplyr::summarise_("mean" = paste("mean(", col_response, ")"),
                    "sd" = paste("sd(", col_response, ")")) 


g.cc.list[[col_response]] <- 
  ggplot(cc.df.sum, 
         aes_string( y = "mean",
                     ymin = "mean - sd",
                     ymax = "mean + sd",
                     x = col_noise)) +
  ggtitle(col_response) +
  ylim(c(0,3)) +
  ylab("Channel capacity") + 
  do.call(theme_jetka, args = plot.args)

if(rep.sample > 5){
  g.cc.list[[col_response]]  <- g.cc.list[[col_response]] + geom_errorbar()
} else {
  g.cc.list[[col_response]]  <- g.cc.list[[col_response]] + geom_point()
}
  


irfmodel.path.list$cc.output.path <- 
  paste(irfmodel.path.list$output.path, cc.id, sep = "/")
dir.create(irfmodel.path.list$cc.output.path, 
           recursive = TRUE)
do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(irfmodel.path.list$cc.output.path, "cc_noise.pdf", sep = "/"),
                           plot = marrangeGrob(grobs = g.cc.list, ncol = 1, nrow = 1))))

saveRDS(
  file = paste(irfmodel.path.list$cc.output.path, "cc_noise.RDS", sep = "/"),
  object = list(
    plots = g.cc.list,
    cc.df = cc.df,
    n.sample  = n.sample,
    rep.sample = rep.sample,
    sd.list = sd.list
  ))   