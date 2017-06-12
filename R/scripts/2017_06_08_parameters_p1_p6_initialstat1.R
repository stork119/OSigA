### ###
### analyse dependence between inital stat1 p1 and p6 parameters
### ###

source("R/optimisation/initialise_optimisation.R")

source("R/model/model_visualisation.R")
source("R/model/model_execution.R")

attach(LoadOptimisationConditions(
  path.optimisation = path.list$optimisation,
  path.optimisation.data = path.list$optimisation.data))

pars <- randomLHS(1000, 3)
pars.i <- 1
pars
parameters.model <- parameters.factor
variables.model <- variables[1:17]
variables.priming.model <-variables.priming[1:17]
no_cores <- 12
registerDoParallel(no_cores)
likelihood.list <- foreach(pars.i = 1:nrow(pars)) %dopar% {
  pars.factor <- 2^(pars[pars.i,]*2 - 1)
  parameters.model[1] <- pars.factor[1]*parameters.factor[c(1)]
  parameters.model[6] <- pars.factor[2]*parameters.factor[c(6)]
  var.factor <- pars.factor[3]
  variables.model[1] <- var.factor*variables[1]
  variables.priming.model[1] <- var.factor*variables.priming[1]
  results <- analyse_model(parameters.model = parameters.model,
                           variables.model  = variables.model,
                           variables.priming.model = variables.priming.model,plot =FALSE, save = FALSE)
  df <- matrix(c(results$likelihood, pars.factor), nrow = 1)
  return(df)
}
stopImplicitCluster()
df.likelihood <- do.call(rbind,likelihood.list)
colnames(df.likelihood) <- c("likelihood", "p1", "p6", "stat1")


pars.factor <- c(1.2665865, 1.3256071, 1.0317962)
parameters.model[1] <- pars.factor[1]*parameters.factor[c(1)]
parameters.model[6] <- pars.factor[2]*parameters.factor[c(6)]
var.factor <- pars.factor[3]
variables.model[1] <- var.factor*variables[1]
variables.priming.model[1] <- var.factor*variables.priming[1]
results <- analyse_model(parameters.model = parameters.model,
                         variables.model  = variables.model,
                         variables.priming.model = variables.priming.model,
                         plot = TRUE,
                         save = TRUE)
write.table(x = df.likelihood %>% data.table() %>% arrange(likelihood),
            file = paste(path.list$optimisation.analysis, "likelihood.csv", sep = "/"), 
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)

df.likelihood <- read.table(file = paste(path.list$optimisation.analysis, "likelihood.csv", sep = "/"), 
                            sep = ",",
                            header = TRUE)
df.likelihood  <- df.likelihood %>% dplyr::mutate(p1p6 =p1/p6, mix = p1*stat1/(p1+p6))



gplot.list[["likelihood_p1p6"]] <- ggplot(data = df.likelihood,
                                          mapping = aes(x = log(p1p6), y = log(stat1), color = likelihood)) + 
  geom_point(size = 2 ) + 
  geom_point(data = df.likelihood[1:100,], aes(x = log(p1p6), y = log(stat1)), color = "red") 
do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.list$optimisation.analysis, "likelihood_p1p6.pdf", sep = "/"),
                           plot = gplot.list[["likelihood_p1p6"]])))
gplot.list[["likelihood_p1p6stat1"]] <- ggplot(data = df.likelihood,
                                               mapping = aes(y = likelihood, x = mix, color = likelihood)) + 
  geom_point(size = 2 ) 
do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.list$optimisation.analysis, "likelihood_p1p6stat1.pdf", sep = "/"),
                           plot = gplot.list[["likelihood_p1p6stat1"]])))

### scatterplot
df.likelihood <- df.likelihood %>% dplyr::mutate


scatterplot_lieklihood(data = df.likelihood,
                       colnames.list = c("p1", "p6", "stat1"),
                       path.list = path.list,
                       filename = "likelihood_scatterplot.pdf")


g <- scatterplot_lieklihood(data = df.likelihood %>% dplyr::mutate(p1_val = parameters.factor[1]*2^p1,
                                                              p6_val = parameters.factor[6]*2^p6,
                                                              stat1_noprm = variables[1]*2^stat1,
                                                              stat1_prm = variables.priming[1]*2^stat1),
                       colnames.list = c("p1_val", "p6_val", "stat1_noprm", "stat1_prm"),
                       path.list = path.list,
                       filename = "likelihood_values_scatterplot.pdf",
                       text_size = 12,
                       width = 16,
                       height = 12)

#### ####
gplot.list[["likelihood_p1p6stat1"]] <- list()
variables[1]
df.likelihood  <- df.likelihood %>% dplyr::mutate(mix = p1*stat1/(stat1+p6),
                                                  mix_var1 = (parameters.factor[1]*2^p1)*(variables[1]*2^stat1)/((variables[1]*2^stat1)+(parameters.factor[6]*2^p6)),
                                                  mix_var2 = (parameters.factor[1]*2^p1)*(variables.priming[1]*2^stat1)/((variables.priming[1]*2^stat1)+(parameters.factor[6]*2^p6)))

gplot.list[["likelihood_p1p6stat1"]][["nopriming_ylim"]] <- 
  ggplot(data = df.likelihood,
         mapping = aes(y = likelihood, x = mix_var1, color = likelihood)) + 
  geom_point(size = 2 )  +
  ylim(c(50000,70000)) +
  xlim(c(0,750)) +
  xlab("p1*stat1/(p6+stat1)") +
  do.call(what = theme_jetka, args = plot.args)  +
  ggtitle("nopriming p1*stat1/(p6+stat1)")

gplot.list[["likelihood_p1p6stat1"]][["priming_ylim"]] <- 
  ggplot(data = df.likelihood,
         mapping = aes(y = likelihood, x = mix_var2, color = likelihood)) + 
  geom_point(size = 2 )  +
  ylim(c(50000,70000)) +
  xlim(c(0,1500)) +
  xlab("p1*stat1/(p6+stat1)") +
  do.call(what = theme_jetka, args = plot.args)  +
  ggtitle("priming p1*stat1/(p6+stat1)")


gplot.list[["likelihood_p1p6stat1"]][["nopriming"]] <- 
  ggplot(data = df.likelihood,
         mapping = aes(y = likelihood, x = mix_var1, color = likelihood)) + 
  geom_point(size = 2 )  +
  xlab("p1*stat1/(p6+stat1)") +
  #ylim(c(50000,70000)) +
  do.call(what = theme_jetka, args = plot.args)  +
  ggtitle("nopriming p1*stat1/(p6+stat1)")

gplot.list[["likelihood_p1p6stat1"]][["priming"]] <- 
  ggplot(data = df.likelihood,
         mapping = aes(y = likelihood, x = mix_var2, color = likelihood)) + 
  geom_point(size = 2 )  +
  xlab("p1*stat1/(p6+stat1)") +
  #ylim(c(50000,70000)) +
  do.call(what = theme_jetka, args = plot.args)  +
  ggtitle("priming p1*stat1/(p6+stat1)")


do.call(what = ggsave,
        args = append(plot.args.ggsave.tmp,
                      list(filename = paste(path.list$optimisation.analysis, "likelihood_p1p6stat1.pdf", sep = "/"),
                           plot = marrangeGrob(grobs = gplot.list[["likelihood_p1p6stat1"]], nrow = 1, ncol = 1))))
### ###
### Verify if absolute values are meaningful
### ###
