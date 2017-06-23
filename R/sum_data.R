# #### ####
# 
# ids <- list.dirs(path.list$optimisation, full.names = FALSE)
# ids <- ids[which(!is.na(as.numeric(ids)))]
# 
# registerDoParallel(no_cores)
# foreach( id = ids[which(!ids %in% optimisation.table$id)] ) %dopar%{
#   try({
#     
#     data.model <-  read.table(
#       file = paste(path.list$optimisation, id, "data_model.csv", sep = "/"),
#       sep = ",",
#       header = TRUE)
# 
#     result <- sapply(fun.likelihood.list,
#                      function(fun.likelihood){
#                        sum( likelihood(
#                          fun.likelihood = fun.likelihood,
#                          data.model = data.model,
#                          data.exp.grouped = data.exp.grouped))
#                      })
# 
#     write.table(x = matrix(result, nrow = 1), paste(path.list$optimisation, id, "optimisation.csv", sep = "/"), sep = ",", row.names = FALSE, col.names = FALSE)
# 
#   })
# }
# stopImplicitCluster()
# 

#### ####
setwd("~/Documents/modelling/")
source("R/optimisation/initialise_optimisation.R")
source("R/sum_data_initialise.R")


gplot.list <- list()
attach(LoadOptimisationConditions(path.optimisation = path.list$optimisation,
                           path.optimisation.data = path.list$optimisation.data))
optimisation.table <- InitiOptimisationTable(path.optimisation = path.list$optimisation,
                                             path.optimisation.data = path.list$optimisation.data)

no_cores <- 16

for( id  in ids[which(!ids  %in% optimisation.table$id)] ){
  try({
        optimisation.table.tmp <-  read_optimisation(path = paste(path.list$optimisation.data, id, sep = "/"),
                                                     id = id, 
                                                     names = names(fun.likelihood.list))
        if(is.null(optimisation.table)){
          optimisation.table <- optimisation.table.tmp
        } else {
          optimisation.table <- rbind(optimisation.table,optimisation.table.tmp)
        }
  })
}

#### ####
fun.likelihood.name <- "sd_data"
optimisation.table <- optimisation.table[order(as.numeric(optimisation.table[,fun.likelihood.name]) ),]
optimisation.best <- c(optimisation.table[order(as.numeric(optimisation.table[,fun.likelihood.name])[1:6]),]$id, "single", "receptors")
#stimulation.list <- (data.exp.grouped %>% ungroup() %>% distinct(stimulation))$stimulation
# data.exp.grouped.all <- data.list$data.exp.norm %>% group_by(priming, stimulation, time) %>% mutate(intensity_sd = var(intensity))
# data.exp.grouped.equal.all <- get_equal_data(data.exp.grouped.all)
# registerDoParallel(no_cores)
# optimisation.table.all <- foreach( i = 1:length(optimisation.best), .combine = rbind) %dopar% {
#   id <- optimisation.best[i]
#   try({
#     parameters <- scan(paste(path.list$optimisation.data, id, "par.txt", sep = "/"))
#     variables <- scan(paste(path.list$optimisation.data, id, "var.txt", sep = "/"))
#     variables.priming <- scan(paste(path.list$optimisation.data, id, "var-priming.txt", sep = "/"))
#     filename <- paste(path.list$optimisation.data, id, "data_model_all_stm.csv", sep = "/")
#     model.simulation <- list()
#     if(file.exists(filename)){
#       model.simulation$data.model <- read.table(
#                   file = filename,
#                   sep = ",",
#                   header = TRUE)
#     } else {
#       model.simulation <- do.call(run_model_mean,
#                                 list(parameters = parameters,
#                                      variables = variables,
#                                      variables.priming = variables.priming,
#                                      tmesh = tmesh,
#                                      tmesh.list = tmesh.list,
#                                      stimulation.list = stimulation.list.all,
#                                      background = background))
#       write.table(x = model.simulation$data.model, 
#                 file = filename,
#                 sep = ",",
#                 row.names = FALSE,
#                 col.names = TRUE)
#     }
#     result <- sapply(fun.likelihood.list,
#                      function(fun.likelihood){
#                        sum( likelihood(
#                          fun.likelihood = fun.likelihood,
#                          data.model = model.simulation$data.model,
#                          data.exp.grouped = data.exp.grouped.optimisation,
#                          data.exp.summarise  = data.exp.summarise.optimisation))
#                      })
#     
#     write.table(x = matrix(result, nrow = 1), 
#                 paste(path.list$optimisation.data, id, "optimisation_all.csv", sep = "/"),
#                 sep = ",", row.names = FALSE, col.names = FALSE)
#     return( matrix(c(id,result), nrow = 1))
#   })
# }
# stopImplicitCluster()

registerDoParallel(no_cores)
data.model.list <- foreach( i = 1:length(optimisation.table$sd_data)) %dopar% {
  id <- optimisation.table$id[i]
  try({
    filename <- paste(path.list$optimisation.data, id, "data_model.csv", sep = "/")
    if(file.exists(filename)){
      data.model <- read.table(
                  file = filename,
                  sep = ",",
                  header = TRUE)
      data.model$type <- id
    return( data.model )
    }
  })
}
stopImplicitCluster()

data.model <- do.call(rbind, data.model.list)

optimisation.table.results <- data.model %>%
  dplyr::group_by(type) %>%
  summarise(likelihood = sum(likelihood)) %>%
  mutate(id = type) %>%
  arrange(likelihood)
#colnames(optimisation.table.all) <- colnames(optimisation.table)
write.table(file = paste(path.list$optimisation.results, "optimisation_ranking_all.csv", sep = ""),
            x = optimisation.table.results,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)

registerDoParallel(no_cores)
data.parameters.list.df <- foreach( i = 1:length(optimisation.table$sd_data)) %dopar% {
  id <- optimisation.table$id[i]
  try({
    filename <- paste(path.list$optimisation.data, id, "par_exp.txt", sep = "/")
    if(file.exists(filename)){
      parameters.id <- scan(file = filename)
      data.parameters <- matrix(parameters.id, nrow = 1)
      colnames(data.parameters) <- paste("p", 1:ncol(data.parameters), sep = "")
      return( data.parameters %>% data.frame() %>% mutate(id = id) )
    }
  })
}
stopImplicitCluster()
data.parameters.df <- do.call(rbind, data.parameters.list.df)

data.parameters.opt <- data.parameters.df[,par.optimised]
#data.parameters.opt$p16 <- data.parameters.opt$p1/data.parameters.opt$p6
data.parameters.opt$id <- data.parameters.df$id

data.parameters.opt.melt <- data.parameters.opt %>% 
  melt(id.vares = id) %>% 
  left_join(optimisation.table.results)
#### save best parameter conditions #####

#parameters[par.optimised] <- parameters.factor[par.optimised]*(parameters.base[par.optimised])^par

par <- as.numeric((data.parameters.df %>% dplyr::filter(id == optimisation.table[1,"id"]))[1,1:10])
par[which(is.na(par))] <- 0
write.table(x = data.frame(factor = parameters.factor*(parameters.base^par),
                           base   = parameters.base,
                           lower  = par.lower,
                           upper  = par.upper),
            file = paste(path.list$optimisation.results, "parameters_conditions.csv", sep = ""),
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")
            
#### ####

gplot.list[["parameters_p1p6"]]  <- ggplot(data = data.parameters %>% 
         dplyr::filter(par %in% c("p1")) %>% 
         left_join((data.parameters %>% dplyr::filter(par %in% c("p6"))), by = "type") %>%
         dplyr::mutate(opt = opt.x/opt.y, id = type) %>%
         left_join(optimisation.table, by = "id"),
       mapping = aes(y = opt, x = sd_data)) + 
  geom_point() + 
  do.call(theme_jetka, args = plot.args)+ 
  ggtitle("p1/p6")
  
do.call(what = ggsave, 
        args = append(plot.args.ggsave, 
                      list(filename = paste(path.list$optimisation.results, "parameters_p1p6.pdf", sep = ""),
                           plot = gplot.list[["parameters_p1p6"]] )))


#### ####
gplot.list[["parameters_optimisation"]] <- ggplot( 
  data.parameters.opt.melt %>% 
    filter(!(id %in% c("single", "receptors"))),
  aes(y = value, x = variable, group = id, colour = likelihood ) ) +
  geom_point() +
  geom_line() + 
  do.call(theme_jetka, plot.args) +
  geom_line(data.parameters.opt.melt %>% filter(likelihood == min(likelihood)), 
            mapping = aes(y = value, x = variable, group = id),
            colour = "red", size = 1.5) + 
  geom_point(data.parameters.opt.melt %>% filter(likelihood == min(likelihood)), 
             mapping = aes(y = value, x = variable, group = id),
             colour = "red", size = 1.5) + 
  ggtitle("Estimated parameters")

do.call(what = ggsave, 
        args = append(plot.args.ggsave, 
                      list(filename = paste(path.list$optimisation.results, "estimated_parameters.pdf", sep = ""),
                           plot = gplot.list[["parameters_optimisation"]])))

#### all estiamtes ####
gplot.list[["optimisation"]] <- list()
  
PlotLocalMin <- function(data,
                         title = ""){
  ggplot(
    data = data %>% 
      filter(!(type %in% c("single", "receptors"))), 
    mapping = aes(x = factor(position, levels = 1:nrow(optimisation.table.results)),
                  y = -log(likelihood), color = likelihood)) + 
    geom_point() + 
    do.call(theme_jetka, plot.args)  +
    ggtitle(title) 
    
}

gplot.list[["optimisation"]][["all"]] <- PlotLocalMin(data = optimisation.table.results %>% 
                                                        dplyr::arrange(likelihood) %>% 
                                                        tibble::rownames_to_column("position"),
                                                      title = "all")
conditions.grid <- expand.grid(
  prm = unique(data.model$priming), 
  stm = unique(data.model$stimulation), 
  tm = unique(data.model$time))

for(conditions.grid.i in 1:nrow(conditions.grid)){
  gplot.list[["optimisation"]][[as.character(conditions.grid.i)]] <-
    data.model %>%
    dplyr::filter(priming == conditions.grid[conditions.grid.i,]$prm,
                  stimulation == conditions.grid[conditions.grid.i,]$stm,
                  time == conditions.grid[conditions.grid.i,]$tm
    ) %>%
    dplyr::arrange(likelihood) %>% 
    tibble::rownames_to_column("position") %>%
    PlotLocalMin(title = paste(conditions.grid[conditions.grid.i,], collapse = "-"))
  
}

do.call(what = ggsave, 
        args = append(plot.args.ggsave, 
                      list(filename = paste(path.list$optimisation.results, "local_minima.pdf", sep = ""),
                           plot = marrangeGrob(grobs = gplot.list[["optimisation"]], ncol = 1, nrow = 1))))
#### parameters table ####
registerDoParallel(no_cores)
data.parameters.list <- foreach( i = 1:length(optimisation.table$sd_data)) %dopar% {
  id <- optimisation.table$id[i]
  try({
    filename <- paste(path.list$optimisation.data, id, "parameters.csv", sep = "/")
    if(file.exists(filename)){
      data.parameters <- read.table(
        file = filename,
        sep = ",",
        header = TRUE)
      data.parameters$par <- paste("p", 1:nrow(data.parameters), sep = "")
      data.parameters$type <- id
      return( data.parameters )
    }
  })
}
stopImplicitCluster()

data.parameters <- do.call(rbind, data.parameters.list) %>% data.table()
data.parameters.opt  <- dcast(data.parameters, type ~ par, value.var = 'opt')
data.parameters.init  <- dcast(data.parameters, type ~ par, value.var = 'init')
data.parameters.init_opt <- left_join(x = data.parameters.init, y = data.parameters.opt, by ="type")
data.parameters.init_opt$id <- data.parameters.init_opt$type
data.parameters.init_opt <- 
  data.parameters.init_opt[,
                           c("id",
                             "type",
                             sapply(1:10,
                                    function(i){
                                      c(
                                        paste("p",i,".","x", sep = ""),
                                        paste("p",i,".","y", sep = ""))}))]

data.parameters.init_opt <- data.parameters.init_opt %>% left_join(optimisation.table, by = "id") %>% dplyr::arrange(sd_data)

write.table(file = paste(path.list$optimisation.results, "parameters.csv", sep = ""),
            x = data.parameters.init_opt,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)

data.parameters.init_opt.melt <- data.parameters.init_opt %>% melt()

params.x <- (data.parameters.init_opt.melt %>% filter(!is.na(value)) %>% distinct(variable))$variable
gplot.list[["parameters"]] <- list()
for(par in 
    unlist(strsplit(x = as.character(params.x[grepl(x = params.x, pattern = "x$")]), split = ".x"))){

  par.compare <- paste(par, c("x", "y"), sep = ".")
  
  gplot.list[["parameters"]][[par]] <- ggplot(data = data.parameters.init_opt.melt %>% 
           dplyr::filter(variable %in% par.compare),
         mapping = aes(x = variable, y = value, group = type )) +
    geom_point() +
    geom_line() + 
    do.call(theme_jetka, plot.args) +
    ggtitle(par)
}

do.call(what = ggsave, 
        args = append(plot.args.ggsave, 
                      list(filename = paste(path.list$optimisation.results, "parameters.pdf", sep = ""),
                           plot = marrangeGrob(grobs = gplot.list[["parameters"]], ncol = 1, nrow = 1))))
#### ####


gplot.list[["models.compare"]] <- list()
for(id in optimisation.table$id[1:10]){
  type.list <- c(id, "single")
  gplot.list[["models.compare"]][[as.character(id)]] <- ggplot(data.model %>% 
         dplyr::filter(type %in% type.list),
       mapping = aes(x = factor(time), y = mean.lmvn, group = type, color = type)) +
  geom_point() +
  geom_line() +
  facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
  geom_point(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
  #geom_line(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
  geom_errorbar(data = data.exp.summarise.optimisation %>% mutate(type = "data"),
                mapping = aes(ymin = mean.lmvn - sqrt(sd.lmvn), 
                              ymax = mean.lmvn + sqrt(sd.lmvn)),
                color = "black") +
  ggtitle(paste("Compare", type.list, collapse = " "))
}
do.call(what = ggsave, 
        args = append(plot.args.ggsave, 
                      list(filename = paste(path.list$optimisation.results, "models_compare.pdf", sep = ""),
                           plot = marrangeGrob(grobs = gplot.list[["models.compare"]], ncol = 1, nrow = 1))))

#### archive ####
# parameters.table.all <- foreach( i = 1:nrow(optimisation.table), .combine = rbind) %do% {
#   id <- optimisation.table$id[i]
#   try({
#     parameters <- scan(paste(path.list$optimisation.data, id, "par.txt", sep = "/"))
#     return( matrix(c(id,parameters), nrow = 1))
#   })
# }
# colnames(parameters.table.all) <- c("id", 1:10)
# write.table(file = paste(path.list$optimisation, "parameters_ranking.csv", sep = ""),
#             x = parameters.table.all,
#             sep = ",",
#             row.names = FALSE,
#             col.names = TRUE)
# write.table(file = paste(path.list$optimisation, "optimisation_ranking.csv", sep = ""),
#             x = optimisation.table,
#             sep = ",",
#             row.names = FALSE,
#             col.names = TRUE)
# 
# plot_results(data = data.exp.grouped.equal.all,
#              path.analysis = path.list$optimisation.data,
#              path.output = path.list$optimisation,
#              data.model.best.list = data.model.list,
#              optimisation.best = optimisation.best[1:5],#st.likelihood[order(st.likelihood$priming),]$id[1:10],
#              plot.title = "all_fullmodel", 
#              filename.data_model = "data_model_all_stm.csv",
#              grid.ncol = 4)
# 
# 
# plot_results(data = data.exp.grouped.equal.all,
#              path.analysis = path.list$optimisation.data,
#              path.output = path.list$optimisation,
#              data.model.best.list = data.model.list,
#              optimisation.best = optimisation.best[1:5],#st.likelihood[order(st.likelihood$priming),]$id[1:10],
#              plot.title = "all", 
#              filename.data_model = "data_model.csv",
#              grid.ncol = 4)
# 
# 
# plot_results(data = data.exp.grouped.equal.all,
#              path.analysis = path.list$optimisation.data,
#              path.output = path.list$optimisation,
#              data.model.best.list = data.model.list,
#              optimisation.best = st.likelihood[order(st.likelihood$priming),]$id[1:10],
#              plot.title = "priming", 
#              filename.data_model = "data_model.csv",
#              grid.ncol = 4)
# 
# 
# 
# 
# #### ####
# # plot_results(data = data.exp.grouped,
# #              path.optimisation = path.list$optimisation,
# #              data.model.best.list = data.model.list,
# #              optimisation.best =  optimisation.table[order(optimisation.table$mvn.mean)[1:5],]$id,
# #              plot.title = "mvn_mean_all", 
# #              filename.data_model = "data_model_all_stm.csv",
# #              grid.ncol = 4)
# # 
# # plot_results(data = data.exp.grouped,
# #              path.optimisation = path.list$optimisation,
# #              data.model.best.list = data.model.list,
# #              optimisation.best = optimisation.table[order(optimisation.table$lmvn)[1:5],]$id,
# #              plot.title = "lmvn")
# # 
# # 
# # parameters.matrix <- matrix(par.def, ncol = 10, nrow = 1)
# # optimisation.best <- optimisation.table[order(optimisation.table$mvn.mean)[1:5],]$id
# # for( i in 1:length(optimisation.best)){
# #   id <- optimisation.best[i]
# #   try({
# #     parameters.matrix <- rbind(parameters.matrix, scan(paste(path.list$optimisation, id, "par.txt", sep = "/")))
# # 
# #   })
# # }
