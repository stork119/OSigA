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

no_cores <- 12

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
optimisation.table <- optimisation.table %>% data.table() %>% dplyr::arrange(sd)
#### ####
fun.likelihood.name <- "sd"

registerDoParallel(no_cores)
data.model.list <- foreach( i = 1:length(optimisation.table$sd)) %dopar% {
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
data.parameters.list.df <- foreach( i = 1:length(optimisation.table$sd)) %dopar% {
  id <- optimisation.table$id[i]
  try({
    filename <- paste(path.list$optimisation.data, id, "parameters_conditions.csv", sep = "/")
    if(file.exists(filename)){
      
      data.parameters <- read.table(file = filename, header = TRUE, sep = ",") %>%
        dplyr::mutate(id = id,
                      value = factor*(base^opt),
                      variable = paste("p", row_number(), sep = "")) %>%
        dplyr::mutate(variable = factor(variable, levels = paste("p", 1:nrow(parameters.conditions), sep = ""))) %>%
        data.table()
      return( data.parameters )
    }
  })
}
stopImplicitCluster()
data.parameters.df <- do.call(rbind, data.parameters.list.df) %>% data.table()

# registerDoParallel(no_cores)
# data.parameters.list.df <- foreach( i = 1:length(optimisation.table$sd)) %dopar% {
#   id <- optimisation.table$id[i]
#   try({
#     filename <- paste(path.list$optimisation.data, id, "par_exp.txt", sep = "/")
#     if(file.exists(filename)){
#       parameters.id <- scan(file = filename)
#       data.parameters <- matrix(parameters.id, nrow = 1)
#       colnames(data.parameters) <- paste("p", par.optimised, sep = "")
#       return( data.parameters %>% data.frame() %>% mutate(id = id) )
#     }
#   })
# }
# stopImplicitCluster()
# data.parameters.df <- do.call(rbind, data.parameters.list.df) %>% data.table()

data.parameters.opt <- data.parameters.df
data.parameters.opt.melt <- data.parameters.opt 

data.parameters.opt.melt <- data.parameters.opt.melt %>% dplyr::left_join(optimisation.table.results)
# data.parameters.opt.melt <- data.parameters.opt %>% 
#   melt(id.vars = "id") %>% 
#   dplyr::left_join(optimisation.table.results) %>%
#   dplyr::filter(variable %in% paste("p", par.optimised, sep = ""))


#### save best parameter conditions #####

#parameters[par.optimised] <- parameters.factor[par.optimised]*(parameters.base[par.optimised])^par

write.table(x = data.parameters.df %>% 
              dplyr::filter(id == optimisation.table[1,"id"]) %>%
              dplyr::mutate(factor = value) %>%
              dplyr::select(factor, base, lower, upper, parameters, parameters.priming),
            file = paste(path.list$optimisation.results, "parameters_conditions.csv", sep = ""),
            col.names = TRUE,
            row.names = FALSE,
            sep = ",")

#### p1p6 ####

# gplot.list[["parameters_p1p6"]]  <- ggplot(data = data.parameters.df %>% 
#          dplyr::filter(par %in% c("p1")) %>% 
#          left_join((data.parameters %>% dplyr::filter(par %in% c("p6"))), by = "type") %>%
#          dplyr::mutate(opt = opt.x/opt.y, id = type) %>%
#          left_join(optimisation.table, by = "id"),
#        mapping = aes(y = opt, x = sd)) + 
#   geom_point() + 
#   do.call(theme_jetka, args = plot.args)+ 
#   ggtitle("p1/p6")
#   
# do.call(what = ggsave, 
#         args = append(plot.args.ggsave, 
#                       list(filename = paste(path.list$optimisation.results, "parameters_p1p6.pdf", sep = ""),
#                            plot = gplot.list[["parameters_p1p6"]] )))
# 
# 
#### estimated parameters how where changed ####
gplot.list[["parameters_optimisation_value"]] <- ggplot( 
  data.parameters.opt.melt %>% 
    dplyr::filter(lower != upper) %>%
    dplyr::filter(!(id %in% c("single", "receptors"))),
  aes(y = log(value), x = variable, group = id, colour = likelihood ) ) +
  geom_point() +
  geom_line() + 
  do.call(theme_jetka, plot.args) +
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
  ggtitle("Estimated parameters log(value)")

do.call(what = ggsave, 
        args = append(plot.args.ggsave, 
                      list(filename = paste(path.list$optimisation.results, "estimated_parameters_value.pdf", sep = ""),
                           plot = gplot.list[["parameters_optimisation_value"]])))

gplot.list[["parameters_optimisation_exp"]] <- ggplot( 
  data.parameters.opt.melt %>% 
    dplyr::filter(lower != upper) %>%
    dplyr::filter(!(id %in% c("single", "receptors"))),
  aes(y = opt, x = variable, group = id, colour = likelihood ) ) +
  geom_point() +
  geom_line() + 
  do.call(theme_jetka, plot.args) +
  geom_line(data.parameters.opt.melt %>%
              dplyr::filter(lower != upper) %>%
              filter(likelihood == min(likelihood)), 
            mapping = aes(y = opt, x = variable, group = id),
            colour = "red", size = 1.5) + 
  geom_point(data.parameters.opt.melt %>% 
               dplyr::filter(lower != upper) %>%
               filter(likelihood == min(likelihood)), 
             mapping = aes(y = opt, x = variable, group = id),
             colour = "red", size = 1.5) + 
  ggtitle("Estimated parameters exponent")

do.call(what = ggsave, 
        args = append(plot.args.ggsave, 
                      list(filename = paste(path.list$optimisation.results, "estimated_parameters_exp.pdf", sep = ""),
                           plot = gplot.list[["parameters_optimisation_exp"]])))


#### all estiamtes ####
gplot.list[["optimisation"]] <- list()
  
PlotLocalMin <- function(data,
                         title = ""){
  ggplot(
    data = data %>% 
      filter(!(type %in% c("single", "receptors"))), 
    mapping = aes(x = factor(position, levels = 1:nrow(optimisation.table.results)),
                  y = log(-likelihood), color = likelihood)) + 
    geom_point() + 
    do.call(theme_jetka, plot.args)  +
    ggtitle(title) 
    
}

gplot.list[["optimisation"]][["all"]] <- PlotLocalMin(data = optimisation.table.results %>% 
                                                        dplyr::arrange(likelihood) %>% 
                                                        dplyr::mutate(position = row_number()),
                                                      title = "all")
conditions.grid <- expand.grid(
  prm = unique(data.model$priming), 
  stm = unique(data.model$stimulation), 
  tm = unique(data.model$time))

for(conditions.grid.i in 1:nrow(conditions.grid)){
  prm <- conditions.grid[conditions.grid.i,]$prm
  stm <- conditions.grid[conditions.grid.i,]$stm
  t <- conditions.grid[conditions.grid.i,]$tm
  gplot.list[["optimisation"]][[as.character(conditions.grid.i)]] <-
    data.model %>%
    dplyr::filter(priming == prm,
                  stimulation == stm,
                  time == t
    ) %>%
    dplyr::arrange(likelihood) %>% 
    dplyr::mutate(position = row_number()) %>%
    PlotLocalMin(title = paste(conditions.grid[conditions.grid.i,], collapse = "-"))
  
}

do.call(what = ggsave, 
        args = append(plot.args.ggsave, 
                      list(filename = paste(path.list$optimisation.results, "local_minima.pdf", sep = ""),
                           plot = marrangeGrob(grobs = gplot.list[["optimisation"]], ncol = 1, nrow = 1))))
#### parameters table ####
# registerDoParallel(no_cores)
# data.parameters.list <- foreach( i = 1:length(optimisation.table$sd)) %dopar% {
#   id <- optimisation.table$id[i]
#   try({
#     filename <- paste(path.list$optimisation.data, id, "parameters.csv", sep = "/")
#     if(file.exists(filename)){
#       data.parameters <- read.table(
#         file = filename,
#         sep = ",",
#         header = TRUE)
#       data.parameters$par <- paste("p", 1:nrow(data.parameters), sep = "")
#       data.parameters$type <- id
#       return( data.parameters )
#     }
#   })
# }
# stopImplicitCluster()
# 
# data.parameters <- do.call(rbind, data.parameters.list) %>% data.table()
# data.parameters.opt  <- dcast(data.parameters, type ~ par, value.var = 'opt')
# data.parameters.init  <- dcast(data.parameters, type ~ par, value.var = 'init')
# data.parameters.init_opt <- left_join(x = data.parameters.init, y = data.parameters.opt, by ="type")
# data.parameters.init_opt$id <- data.parameters.init_opt$type
# data.parameters.init_opt <- 
#   data.parameters.init_opt[,
#                            c("id",
#                              "type",
#                              sapply(1:length(parameters.factor),
#                                     function(i){
#                                       c(
#                                         paste("p",i,".","x", sep = ""),
#                                         paste("p",i,".","y", sep = ""))}))]
# 
# data.parameters.init_opt <- data.parameters.init_opt %>% left_join(optimisation.table, by = "id") %>% dplyr::arrange(sd)
# 
# write.table(file = paste(path.list$optimisation.results, "parameters.csv", sep = ""),
#             x = data.parameters.init_opt %>% dplyr::arrange(sd),
#             sep = ",",
#             row.names = FALSE,
#             col.names = TRUE)
# 
#   data.parameters.init_opt.melt <- data.parameters.init_opt %>% reshape2::melt(id.vars = c("id", "type", "sd_data", "sd", "data")) %>% data.table()
# 
# params.x <- (data.parameters.init_opt.melt %>% filter(!is.na(value)) %>% distinct(variable))$variable
# gplot.list[["parameters"]] <- list()
# for(par in 
#     unlist(strsplit(x = as.character(params.x[grepl(x = params.x, pattern = "x$")]), split = ".x"))){
# 
#   par.compare <- paste(par, c("x", "y"), sep = ".")
#   
#   gplot.list[["parameters"]][[par]] <- ggplot(data = data.parameters.init_opt.melt %>% 
#            dplyr::filter(variable %in% par.compare),
#          mapping = aes(x = variable, y = value, group = type )) +
#     geom_point() +
#     geom_line() + 
#     do.call(theme_jetka, plot.args) +
#     ggtitle(par)
# }
# 
# do.call(what = ggsave, 
#         args = append(plot.args.ggsave, 
#                       list(filename = paste(path.list$optimisation.results, "parameters.pdf", sep = ""),
#                            plot = marrangeGrob(grobs = gplot.list[["parameters"]], ncol = 1, nrow = 1))))


write.table(file = paste(path.list$optimisation.results, "parameters_ranking.csv", sep = ""),
            x = data.parameters.df %>% 
              reshape2::dcast(id ~ variable, value.var = "value") %>% 
              left_join(y = optimisation.table, by = "id") %>%
              dplyr::arrange(sd),
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)

write.table(file = paste(path.list$optimisation.results, "parameters_ranking_opt.csv", sep = ""),
            x = data.parameters.df %>% 
              reshape2::dcast(id ~ variable, value.var = "opt") %>% 
              left_join(y = optimisation.table, by = "id") %>%
              dplyr::arrange(sd),
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)

#### ####


gplot.list[["models.compare_log"]] <- list()
for(id in optimisation.table$id[1:20]){
  type.list <- c(id, "single")
  gplot.list[["models.compare_log"]][[as.character(id)]] <- ggplot(data.model %>% 
         dplyr::filter(type %in% type.list),
       mapping = aes(x = factor(time), 
                     y = mean.lmvn,
                     group = type,
                     color = type,
                     ymin = mean.lmvn - sqrt(sd.lmvn), 
                     ymax = mean.lmvn + sqrt(sd.lmvn))) +
  geom_point() +
  geom_line() +
    geom_errorbar() +
  facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
  geom_point(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
  #geom_line(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
  geom_errorbar(data = data.exp.summarise.optimisation %>% mutate(type = "data"),
                color = "black") +
  ggtitle(paste("Compare", type.list, collapse = " "))
}
do.call(what = ggsave, 
        args = append(plot.args.ggsave, 
                      list(filename = paste(path.list$optimisation.results, "models_compare_log.pdf", sep = ""),
                           plot = marrangeGrob(grobs = gplot.list[["models.compare_log"]], ncol = 1, nrow = 1))))


gplot.list[["models.compare"]] <- list()
for(id in optimisation.table$id[1:20]){
  type.list <- c(id, "single")
  gplot.list[["models.compare"]][[as.character(id)]] <- ggplot(data.model %>% 
                                                                 dplyr::filter(type %in% type.list),
                                                               mapping = aes(x = factor(time), 
                                                                             y = m.norm,
                                                                             group = type,
                                                                             color = type,
                                                                             ymin = m.norm - sqrt(sd.norm), 
                                                                             ymax = m.norm + sqrt(sd.norm))) +
    geom_point() +
    geom_line() +
    geom_errorbar() +
    facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
    geom_point(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
    #geom_line(data = data.exp.summarise.optimisation %>% mutate(type = "data"), color = "black") +
    geom_errorbar(data = data.exp.summarise.optimisation %>% mutate(type = "data"),
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
