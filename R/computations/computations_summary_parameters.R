### ###
### script parameters correlation
### ###

no_cores = 12
optimisation.table <- read.table(paste(path.list$optimisation.results, "otimisation_fit.csv", sep = "/" ),
                                 header = TRUE,
                                 sep = ",")
optimisation.initiate <- InitiateOptimisation(
  path.list = path.list)
parameters.conditions <- optimisation.initiate$parameters.conditions
parameters.conditions$variable <- paste("p",
                                        1:nrow(parameters.conditions),
                                        sep ="")
registerDoParallel(no_cores)
parameters.estimated.list <- 
  foreach( i = 1:length(optimisation.table$par.id)) %dopar% {
    par_id <-  optimisation.table[i,]$par.id
    parameters.tmp <- read.table(
      file = paste(path.list$optimisation.data, par_id, "parameters.csv", sep = "/"),
      header =  TRUE,
      sep = ","
    )
    parameters.tmp <- parameters.tmp %>%
      dplyr::mutate(data.id = id) %>%
      dplyr::mutate(par.id = par_id) %>%
      dplyr::select(-id)
  }
stopImplicitCluster()
parameters.estimated <- 
  do.call(what = rbind,
          args = parameters.estimated.list)

#install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")

parameters.estimated <- parameters.estimated %>% 
  dplyr::filter(likelihood < -50000)

parameters.estimated <- parameters.estimated %>% 
  #dplyr::mutate(p1p6 = p11*p21*p1/(p6+p19*p11))
  dplyr::mutate(p1p6 = p1/p6)

df <- parameters.estimated[, c(paste("p",c(2,3,8,9,10), sep = ""), "p1p6")]

M <- cor(df)
#M <- M[-c(4,5,7),-c(4,5,7)]
do.call(pdf, 
        args = append(plot.args.ggsave,
                      list(file = 
                             paste(path.list$optimisation.results,
                                   "parameters_correlations.pdf", sep ="/"))))
chart.Correlation(log(df))
corrplot(M, method = "color")
corrplot(M, method = "number")
dev.off()

parameters.estimated <- parameters.estimated %>% 
  dplyr::mutate(p1p6 = p1/p6) %>%
  dplyr::mutate(p9p10 = p9/p10) %>%
  dplyr::mutate(p8p9 = p8*p9*p10) 
parameters.estimated.melt <- parameters.estimated %>% 
  melt(id.vars = c("likelihood", "par.id", "data.id")) %>%
  data.table()
parameters.estimated.melt <- parameters.estimated.melt %>% dplyr::mutate(variable = as.character(variable))

par <- "p1p6"
q <- quantile(parameters.estimated.melt$likelihood, probs = c(0.75))[[1]]
colnames(df)
gplot.list <- list()
for(par in colnames(df)){
gplot.list[[par]] <- ggplot(
  data = parameters.estimated.melt %>%
    dplyr::filter(variable %in% c(par),
                  likelihood < q),
  aes(
    x = log(value),
    y = log(-likelihood))
) +
  geom_point() +
  #xlim(c(-16, -15)) + 
  xlab(paste("log(",par,")")) +
  ggtitle(paste("Dependence of likelihood on parameters (",par,")"))
}


do.call(what = ggsave,
                  args = append(plot.args.ggsave,
                                list(filename = paste(path.list$optimisation.results,
                                                      "parameters_dependence.pdf", sep ="/"),
                                     plot = marrangeGrob(grobs = gplot.list, ncol = 1, nrow = 1))))

