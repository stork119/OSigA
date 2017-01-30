# setwd("~/Documents/modelling/")
# source("R/libraries.R")
# source("R/initialise.R")

#### ####
read_optimisation <- function(path, id){
  optimisation <- read.table(file = paste(path, "optimisation.csv", sep = "/"),
                             sep = ",", header = FALSE)
  data.table(
    id = id,
    mvn.mean = optimisation[1,1],
    mvn = optimisation[1,2],
    lmvn = optimisation[1,3]
  )
}


plot_results <- function(data = data.exp.grouped,
                         path.analysis,
                         data.model.best.list = data.model.list,
                         optimisation.best,
                         grid.ncol = 2,
                         grid.nrow = 2,
                         plot.title = "opt",
                         filename.data_model = "data_model.csv"){
  plot.list <- list()
  
  for( i in 1:length(optimisation.best) ){
    id <- optimisation.best[i]
    try({
      data.model.best.list[[as.character(id)]] <- read.table(
        file = paste(path.analysis, id, filename.data_model, sep = "/"),
        sep = ",",
        header = TRUE)
      plot.title.id <- paste(plot.title, "place:", i, ", id:", id, ",", sep = "")
      plot.list[[as.character(id)]] <- compare_models(
        data = data,
        data.model.list = data.model.best.list[c("single", "double", as.character(id))],
        plot.title = paste(plot.title, "place:", i, ", id:", id, ",", sep = ""),
        filename = paste(path.analysis, plot.title, "_", i, ".pdf", sep = ""),
        plot.save = FALSE,
        plot.return = TRUE
      )
      
    })
  }
  pdf(file = paste(path.analysis, plot.title, ".pdf", sep = ""),
      width = 20,
      height = 12,
      useDingbats = FALSE)
  print(marrangeGrob(unlist(plot.list, recursive = FALSE), ncol = grid.ncol, nrow =  grid.nrow))
  dev.off()
}
#### ####
path.analysis <- aste(path.output, "cmaes/mvn/2017-01-28/", sep = "/")

optimisation.table <- data.table(
  id = character(),
  mvn.mean = numeric(),
  mvn = numeric(),
  lmvn = numeric()
)

path.single <- paste(path.optimisation, "single", sep = "/")
optimisation.table <- rbind(optimisation.table, read_optimisation(path = path.single, id = "single"))

path.receptors <- paste(path.optimisation, "receptors", sep = "/")
optimisation.table <- rbind(optimisation.table, read_optimisation(path = path.receptors, id = "receptors"))

