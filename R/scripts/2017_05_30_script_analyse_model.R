### ###
### 2017-05-30-script-run_model
### ###

source("R/initialise.R")
source("R/graphics/libraries.R")
# 2017_05_17 initialise
#### read parameters from command line ####
parameters <- scan()
write.table(x = parameters,
            file = "/media/knt/D/KN/ClustSense/Modelling/lib/StochSense/StochSens_Ver2.3/models/JAKSTAT/Input/parameters.csv", 
            sep =",", 
            col.names =FALSE, 
            row.names = FALSE)
write.table(x = paste("p", 1:length(parameters), sep = ""),
            file = "/media/knt/D/KN/ClustSense/Modelling/lib/StochSense/StochSens_Ver2.3/models/JAKSTAT/Input/parameters-names.csv", 
            sep =",", 
            col.names =FALSE, 
            row.names = FALSE)
write.table(x = variables[1:17],
            file = "/media/knt/D/KN/ClustSense/Modelling/lib/StochSense/StochSens_Ver2.3/models/JAKSTAT/Input/variables.csv", 
            sep =",", 
            col.names =FALSE, 
            row.names = FALSE)

write.table(x = variables.priming[1:17],
            file = "/media/knt/D/KN/ClustSense/Modelling/lib/StochSense/StochSens_Ver2.3/models/JAKSTAT/Input/variables-priming.csv", 
            sep =",", 
            col.names =FALSE, 
            row.names = FALSE)
write.table(x = paste("y", 1:length(parameters), sep = ""),
            file = "/media/knt/D/KN/ClustSense/Modelling/lib/StochSense/StochSens_Ver2.3/models/JAKSTAT/Input/variables-names.csv", 
            sep =",", 
            col.names =FALSE, 
            row.names = FALSE)

#### ####
data.model.list[["model_lmvn"]] <- 
  run_model_mean(parameters = parameters, 
               variables  = variables,
               variables.priming = variables.priming,
               tmesh = tmesh,
               tmesh.list = tmesh.list,
               stimulation.list = stimulation.list,
               background = background,
               time_interval = 100,
               time_computation = 1000*60*5)




g.list[["models_comaprison"]] <- 
  do.call(what = plot_points,
          args = append(plot.args,
                        list(
                          data = data.list$data.exp.summarise %>% 
                            dplyr::mutate(yerror = sqrt(sd.lmvn)),
                          x = "time",
                          y = "mean.lmvn",
                          yerror = "yerror",
                          facet_grid_group_y = "stimulation",
                          facet_grid_group_x = "priming",
                          ylim_min = NULL,
                          ylim_max = NULL
                        ))) + 
  geom_line(data = data.model.list[["model_lmvn"]]$data.model, 
            mapping = 
              aes(x = factor(time),
                  y = log(m.norm),
                  group = interaction(priming, stimulation)),
            color = "red") +
  geom_line(data = data.model.list$single, 
            mapping = 
              aes(x = factor(time),
                  y = log(m.norm),
                  group = interaction(priming, stimulation)),
            color = "green") +
  geom_line(data = data.model.list$double, 
            mapping = 
              aes(x = factor(time),
                  y = log(m.norm),
                  group = interaction(priming, stimulation)),
            color = "blue")
ggsave(plot = g.list[["models_comaprison"]],
       filename = paste(path.data.output, 
                        "models_comaprison.pdf", 
                        sep = "/" ),
       width  = plot.args$width, 
       height = plot.args$height, 
       useDingbats = plot.args$useDingbats)  
#### likelihood ####

data <- 
  data.model.list[["model_lmvn"]]$data.model %>% 
    left_join(data.list$data.exp.summarise, by =  c("time", "priming", "stimulation")) %>%
    dplyr::mutate(likelihood = ComputeLikelihood.lmvn.nusigma(nu = mean.lmvn.y, sigma = sd.lmvn.y, X = m.norm.x)) %>% 
  data.table()
  
fun.likelihood <- ComputeLikelihood.lmvn.nusigma

g.list[["model_likelihood"]] <- list()
PlotModelLikelihood <- function(fun.likelihood){
  do.call(what = plot_points,
          args = append(plot.args,
                        list(
                          data = data.model.list[["model_lmvn"]]$data.model %>% 
                            left_join(data.list$data.exp.summarise, by =  c("time", "priming", "stimulation")) %>%
                            dplyr::mutate(likelihood = fun.likelihood(nu = mean.lmvn.y, sigma = sd.lmvn.y, X = m.norm.x)) %>% 
                            data.table(),
                          x = "time",
                          y = "likelihood",
                          facet_grid_group_y = "stimulation",
                          facet_grid_group_x = "priming",
                          ylim_min = NULL,
                          ylim_max = NULL
                        ))) +
  geom_point(data = data.model.list[["single"]] %>% 
               left_join(data.list$data.exp.summarise, by =  c("time", "priming", "stimulation")) %>%
               dplyr::mutate(likelihood = fun.likelihood(nu = mean.lmvn.y, sigma = sd.lmvn.y, X = m.norm.x)) %>% 
               data.table(), 
             mapping = aes(x = factor(time), y = likelihood), color ="green") +
  geom_point(data = data.model.list[["double"]] %>% 
               left_join(data.list$data.exp.summarise, by =  c("time", "priming", "stimulation")) %>%
               dplyr::mutate(likelihood = fun.likelihood(nu = mean.lmvn.y, sigma = sd.lmvn.y, X = m.norm.x)) %>% 
               data.table(), 
             mapping = aes(x = factor(time), y = likelihood), color ="blue")
}

g.list[["model_likelihood"]][["lmvn"]] <- PlotModelLikelihood(fun.likelihood = ComputeLikelihood.lmvn.nusigma) + ggtitle("lmvn")
g.list[["model_likelihood"]][["bias"]] <- PlotModelLikelihood(fun.likelihood = ComputeLikelihood.lmvn.bias) + ggtitle("bias")
g.list[["model_likelihood"]][["rse"]] <- PlotModelLikelihood(fun.likelihood = ComputeLikelihood.lmvn.rse) + ggtitle("rse")

ggsave(plot = marrangeGrob(grobs = g.list[["model_likelihood"]], ncol = 1, nrow =1),
       filename = paste(path.data.output, 
                        "model_likelihood.pdf", 
                        sep = "/" ),
       width  = plot.args$width, 
       height = plot.args$height, 
       useDingbats = plot.args$useDingbats)
#### optimisation ####

data.model.list[["model_lmvn"]]$data.model$likelihood  <- 
  likelihood(data.model = data.model.list[["model_lmvn"]]$data.model,
             data.exp.grouped = data.list$data.exp.norm,
             data.exp.summarise =  data.list$data.exp.summarise,
             fun.likelihood = fun.likelihood.list$sd_data)
                                               

data.model.list$single$likelihood  <- 
  likelihood(data.model = data.model.list$single,
             data.exp.grouped = data.list$data.exp.norm,
             data.exp.summarise =  data.list$data.exp.summarise,
             fun.likelihood = fun.likelihood.list$sd_data)


data.model.list$double$likelihood  <- 
  likelihood(data.model = data.model.list$double,
             data.exp.grouped = data.list$data.exp.norm,
             data.exp.summarise =  data.list$data.exp.summarise,
             fun.likelihood = fun.likelihood.list$sd_data)

sum(data.model.list[["model_lmvn"]]$data.model$likelihood )
sum(data.model.list$single$likelihood)
sum(data.model.list$double$likelihood)

data.model.list[["model_lmvn"]]$data.model$type <- "model_lmvn"
data.model.list$single$type <- "single"
data.model.list$double$type <- "double"

data <- rbind(data.model.list[["model_lmvn"]]$data.model, data.model.list$single)
data <- rbind(data, data.model.list$double)


ggplot(data = data %>% dplyr::filter(stimulation != 5),
       mapping = aes(x = factor(time), y = likelihood, color = type)) +
  geom_point() + 
  do.call(theme_jetka, args = plot.args) + 
  facet_grid(priming ~ stimulation)
#### ####
saveRDS(object = 
          list(
            data.model.list = data.model.list, 
            data.list = data.list,
            g.list = g.list),
        file = paste(path.output,
                     "rds",
                     "2017-06-01.rds",
                     sep = "/"))