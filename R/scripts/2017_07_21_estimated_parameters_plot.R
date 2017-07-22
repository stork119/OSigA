

setwd("~/Documents/modelling/")
source("R/optimisation/initialise_optimisation.R")
source("R/sum_data_initialise.R")
source("R/model/ut/model_ut.R")

gplot.list <- list()
optimisation.conditions.toload <-
  LoadOptimisationConditions(path.optimisation = path.list$optimisation,
                             path.optimisation.data = path.list$optimisation.data,
                             maxit.tmp = Inf)
rm(list = labels(optimisation.conditions.toload))
attach(optimisation.conditions.toload)

#### estimated parameters 

par.optimised
plots.parameters <- list()
plots.parameters$nonpriming <- which(which(parameters.conditions$parameters != 0) %in% par.optimised)
plots.parameters$priming <- which(which(parameters.conditions$parameters.priming != 0) %in% par.optimised)


data.parameters.opt.melt.2 <-
  data.parameters.opt.melt %>%
  dplyr::filter(likelihood < 0) %>%
  dplyr::mutate(gradient = log(-likelihood))

g <- quantile(data.parameters.opt.melt.2$gradient, probs = 0.1)

data.parameters.opt.melt.2 <-
  data.parameters.opt.melt.2 %>%
  dplyr::filter(gradient > g)
  
data.lower <- parameters.conditions %>%
  mutate(
    id = paste("p", row_number(), sep = ""),
    variable = paste("p", row_number(), sep = ""),
    id = "lower",
    value = factor*base^lower,
    gradient = 0,
    likelihood = 0)

data.upper <- parameters.conditions %>%
  mutate(
    variable = paste("p", row_number(), sep = ""),
    id = "upper",
    value = factor*base^upper,
    gradient = 0,
    likelihood = 0)

ggplot( 
  data.parameters.opt.melt.2 %>%
    dplyr::filter(likelihood < 0) %>%
    dplyr::filter(lower != upper) %>%
    dplyr::filter(!(id %in% c("single", "receptors"))),
  aes(y = log(value), 
      x = variable, 
      group = id, 
      colour = gradient ) ) +
  geom_point() +
  geom_line() + 
  do.call(theme_jetka, plot.args) +
  
  scale_colour_gradient(low = "white", high = "black") +
  geom_line(data.lower   %>%
              dplyr::filter(lower != upper),
            mapping = aes(y = log(value), 
                                      x = variable),
            colour = "yellow", size = 1.5) +
  geom_line(data.upper   %>%
              dplyr::filter(lower != upper),
            mapping = aes(y = log(value), 
                          x = variable),
            colour = "yellow", size = 1.5) +
  geom_line(data.parameters.opt.melt %>%
              dplyr::filter(lower != upper) %>%
              filter(likelihood == min(likelihood)), 
            mapping = aes(y = log(value), x = variable, group = id),
            colour = "red", size = 1.5) + 
  geom_point(data.parameters.opt.melt %>% 
               dplyr::filter(lower != upper) %>%
               filter(likelihood == min(likelihood)), 
             mapping = aes(y = log(value), x = variable, group = id),
             colour = "red", size = 1.5) + 
  scale_colour_gradient(low = "white", high = "black") +
  ggtitle("Estimated parameters log(value)")


log(-data.parameters.opt.melt$likelihood)