### ###
###
### ###


source(file = "R/optimisation/initialise_optimisation.R")

likelihood.cv.list <- list()
# 
# path.list <-
#   LoadOptimisationPaths(
#     path.output = "resources/output/",
#     id = "2017-10-17-1"
#   )
# 
# likelihood.cv.list[["0"]] <-
#   read.table(paste(path.list$optimisation.results, "otimisation_cv.csv", sep = "/"), 
#              sep = ",",
#              header = TRUE)
# likelihood.cv.list[["0"]]$par <- 0


path.list <-
  LoadOptimisationPaths(
    path.output = "resources/output/",
    id = "2017-10-17-2"
  )
for( par in parameters.list[c(1,2,3,4)]){
  #par <- 2
  path.list.par <-
    LoadOptimisationPaths(
      path.output = paste("resources/output/", sep = "/"),
      id.list = list(path.list$id, par)
    )
    

  likelihood.cv.list[[as.character(par)]] <-
    read.table(paste(path.list.par$optimisation.results, "otimisation_cv.csv", sep = "/"), 
             sep = ",",
             header = TRUE)[1:10,]
  likelihood.cv.list[[as.character(par)]]$par <- par
  
  
}

#saveRDS(likelihood.cv.list, file = "resources/output/optimisation/2017-10-17-2/analysis/likelihood.RDS")

  likelihood.cv.df <- do.call(what = rbind, args = likelihood.cv.list)  
  likelihood.cv.df <- likelihood.cv.df[-2,]
  likelihood.cv.df <- likelihood.cv.df %>% dplyr::filter(likelihood < -51000)
  ggplot(data = likelihood.cv.df, 
         aes(x = log(-likelihood), 
             color = factor(par), 
             group = factor(par))) + 
    geom_density() + 
    xlim(c(10.75, 11))
  
  
  
  ks.0 <- (likelihood.cv.df %>% dplyr::filter(par == 0))$likelihood
  ks.1 <- (likelihood.cv.df %>% dplyr::filter(par == 4))$likelihood
  ks.test(ks.0, ks.1)
  
  
  ks.0.df <- data.frame(id = (0:(length(ks.0)-1))/length(ks.0), likelihood = sort(ks.0), par = 0)
  ks.1.df <- data.frame(id = (0:(length(ks.1)-1))/length(ks.1), likelihood = sort(ks.1), par = 1)
  ggplot(rbind(ks.0.df, ks.1.df), mapping = aes(x = id, y = log(-likelihood), group = par, color = par)) + geom_line()
# }