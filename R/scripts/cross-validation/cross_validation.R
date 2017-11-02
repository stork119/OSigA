### ###
### cross validation ###
### ###


### ###
### 2017-07-11-sigmapoints-script
### ###

#setwd("~/Documents/modelling/")
#### libraries ####
source("R/optimisation/initialise_optimisation.R")
source("R/sum_data_initialise.R")
source("R/model/ut/model_ut.R")
source("R/scripts/2017_07_11_sigmapoints_library.R")
source("R/jakstat_model.R")
#### computations ####
model.list <- c("modelfull", "model_p1", "model_p10")

model.parameters.dt <- data.table(
  model     = c("modelfull", "model_p1", "model_p10"),
  theta_num = c(16, 17, 17)
)

analyse_name = "cross-validation-3"

fun.likelihood <- fun.likelihood.list.sd
data_sample_size.list <- 250#seq(from = 250, to = 2000, by = 250)
cv_times <- 10000
no_cores <- 4

model <- model.list[1]
likelihood.df.list <- list()
#### cross_validate hypotesis ####
for(model in model.list[c(1,2,3)]){
  
  path.list <- LoadOptimisationPaths(
    path.output = "resources/output/",
    id = paste("2017-07-19-ut-", model, sep = "")
  )
  
  path.optimisation <- path.list$optimisation
  path.optimisation.data <- path.list$optimisation.data
  path.optimisation.analysis <- path.list$optimisation.analysis
  path.optimisation.results <- path.list$optimisation.results
  
  gplot.list <- list()
  optimisation.conditions.toload <-
    LoadOptimisationConditions(optimisation.path =  path.list$optimisation)
  rm(list = labels(optimisation.conditions.toload))
  attach(optimisation.conditions.toload)
  
  
  sigmapoints <- LoadSigmapointsConditions(path.optimisation = path.list$optimisation)
  #data.model <- read.table(paste("resources/output//optimisation/2017-07-19-ut-model_p1///data/", results.id, "data_model.csv", sep = "/"), header = TRUE, sep = ",")
  
  variables.model <- scan(paste(path.list$optimisation, "variables.csv", sep = "/"))
  variables.priming.model <- scan(paste(path.list$optimisation, "variables-priming.csv", sep = "/"))
  ##
  #results.id <- "108"
  # parameters.df <- read.table(paste(path.list$optimisation, "parameters_conditions.csv", sep = "/"), header = TRUE, sep = ",") # powinno być results
  # parameters.df <- parameters.df %>% 
  #   # dplyr::mutate(par = opt) %>%
  #   # dplyr::select(-c(opt, init)) %>%
  #   # dplyr::mutate(factor = factor*base^par) %>% 
  #   dplyr::mutate(par = 0)
  
  ## no results 
  parameters.df <- parameters.conditions
  parameters.df <- parameters.df %>% 
    dplyr::mutate(par = 0)
  
  
  res <- analyse_model_ut(parameters = parameters.factor,
                          variables.model = variables.model,
                          variables.priming.model = variables.priming.model,
                          sigmapoints = sigmapoints,
                          analyse_name = analyse_name,
                          model.computations = list(raw = TRUE, priming = TRUE),
                          parameters.df = parameters.df,
                          fun_modify_parameters = PrepareModelParameters.ut,
                          fun_modify_input = PrepareModelArguments.ut.multiple,
                          parameters.conditions = parameters.conditions,
                          fun_parameters_penalty = fun_parameters_penalty_sigmapoints
  )
  compare_distribution(data.exp.grouped.optimisation = data.exp.grouped.optimisation, data.model = res$data.model, analyse_name = res$analyse_name)
  
  path <- res$path
  data.model <- res$data.model
  
  likelihood.list <- list()
  likelihood.list.add <- list()
  for(sample_size in data_sample_size.list){
    print(sample_size)
    registerDoParallel(no_cores)
    likelihood.list.tmp <- 
      foreach(i = 1:cv_times) %dopar% {
        
        data.cv.list <- list()
        data.cv.list$grouped <- get_equal_data(data = data.list$data.exp, 
                                               sample_size = sample_size)
        
        data.cv.list$summarise  <- 
          data.cv.list$grouped %>% 
          dplyr::group_by(priming,
                          stimulation,
                          time) %>%
          dplyr::summarise(m.norm = mean(intensity),
                           sd.norm   = var(intensity)) %>%
          dplyr::mutate(mean.lmvn = lmvn.mean(m = m.norm, sd = sd.norm),
                        sd.lmvn = lmvn.sd(m = m.norm, sd = sd.norm))
        likelihood.add <- sample_size*log(2*pi) + sum(data.cv.list$grouped$logintensity)
        likelihood.sd <- sum(likelihood(data.model = data.model,# %>% dplyr::filter(time %in% tmesh[tmesh.list]),
                                        data.exp.grouped = data.cv.list$grouped,
                                        data.exp.summarise =   data.cv.list$summarise,
                                        fun.likelihood = fun.likelihood))
        return(list(sd = likelihood.sd, add= likelihood.add))
      }
    stopImplicitCluster()
    likelihood.list[[as.character(sample_size)]] <- sapply(likelihood.list.tmp, function(l){l$sd})
    likelihood.list.add[[as.character(sample_size)]] <- sapply(likelihood.list.tmp, function(l){l$add})
  }
  
  likelihood.df <- matrix(nrow = length(data_sample_size.list)*cv_times, ncol = 4)
  for(i in 1:length(data_sample_size.list)){
    likelihood.df[((i-1)*cv_times + 1):(i*cv_times),1] <- model
    likelihood.df[((i-1)*cv_times + 1):(i*cv_times),2] <- data_sample_size.list[i]
    likelihood.df[((i-1)*cv_times + 1):(i*cv_times),3] <- likelihood.list[[i]]
    likelihood.df[((i-1)*cv_times + 1):(i*cv_times),4] <- likelihood.list.add[[i]]
    
  }
  colnames(likelihood.df) <- c("model", "sample_size", "likelihood", "likelihood.add")
  likelihood.df <- data.table(likelihood.df) %>%
    mutate(likelihood = as.numeric(likelihood),
           likelihood.add = as.numeric(likelihood.add),
           sample_size = as.numeric(sample_size))
  
  
  write.table(x = likelihood.df,
              file = paste(path, "cross-validation.csv", sep = "/"),
              sep = ",",
              row.names = FALSE,
              col.names = TRUE)
  likelihood.df.list[[model]] <- likelihood.df
}


likelihood.df.all <- do.call(rbind, likelihood.df.list)

#### plotting ####
results <- saveRDS(object = likelihood.df.all, file = "resources/output/cross_validation/cross-validation-3/cross_valid.RDS")

results <- readRDS("resources/output/cross_validation/cross-validation-3/cross_valid.RDS")

likelihood.df.all %>% 
  left_join(model.parameters.dt, by = "model") %>%
  dplyr::mutate(aic = as.numeric(likelihood.add) + likelihood  + 2*theta_num,
                bic = as.numeric(likelihood.add) + likelihood  +
                  log(sample_size*84)*theta_num/2)

ggplot(likelihood.df.all %>% 
         dplyr::group_by(model, sample_size) %>%
         dplyr::summarise(mean = mean(-likelihood),
                          sd = var(-likelihood),
                          logmean = mean(log(-likelihood)),
                          logsd = var(log(-likelihood))),
       mapping =
         aes(x = sample_size,
             y = mean,
             ymin = mean - sqrt(sd),
             ymax = mean + sqrt(sd),
             group = model,
             color = model,
             fill  = model)) +
  geom_line(alpha = 0.5) + 
  geom_ribbon(alpha = 0.5) +
  do.call(theme_jetka, args = plot.args)

likelihood.df.all.2 <- likelihood.df.all
likelihood.df.all.2$model <- "model_all"
likelihood.df.all <- likelihood.df.all %>%
  rbind(likelihood.df.all.2)
models.name.list <- 
  list(
    "modelfull" = 0,
    "model_p1" = 1, 
    "model_p4p10" = 4, 
    "model_p10" = 10,
    "model_all" = 2)

likelihood.df.all <-
  likelihood.df.all %>%
  dplyr::mutate(model = as.character(model)) %>%
  dplyr::mutate(model_new = models.name.list[as.character(model)][[1]])

g <- ggplot(likelihood.df.all %>% dplyr::filter(sample_size == 250, model != "model_p4"),
       aes(x = likelihood,#log(-likelihood), 
           color = factor(model), 
         #  fill = model,
           alpha = 0.0001)) +
  geom_density() +
  theme_jetka()
do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = "resources/output/cross_validation/cross-validation-3/densities-poster-2017-10-24.pdf",
                           plot = g)))


  quantile(ks.1, c(.99, .1))
###
ks.1 <- (likelihood.df.all %>% dplyr::filter(sample_size == 250, model == "modelfull"))$likelihood
ks.2 <- (likelihood.df.all %>% dplyr::filter(sample_size == 250, model == "model_p1"))$likelihood 
ks.2 <- (likelihood.df.all %>% dplyr::filter(sample_size == 250, model == "model_p4"))$likelihood 
ks.2 <- (likelihood.df.all %>% dplyr::filter(sample_size == 250, model == "model_p10"))$likelihood 
ks.2 <- (likelihood.df.all %>% dplyr::filter(sample_size == 250, model == "model_p4p10"))$likelihood
ks.2 <- (likelihood.df.all %>% dplyr::filter(sample_size == 250, model == "model_all"))$likelihood
ks.test(x = ks.1, y = ks.2)
ks.df <- data.frame(ks1 = sort(ks.1), ks2 = sort(ks.2))
ks.df$x <- 1:nrow(ks.df)
ggplot(ks.df, aes(x= x, y = ks1)) + geom_point() + geom_point(aes(x= x, y = ks2))