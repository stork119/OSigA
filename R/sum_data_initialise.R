# setwd("~/Documents/modelling/")
# source("R/libraries.R")
# source("R/initialise.R")

#### ####
optimisation.table <- data.table(
  id = character(),
  mvn.mean = numeric(),
  mvn = numeric(),
  lmvn = numeric(),
  mvn.sd_const = numeric(),
  lmvn.data = numeric(),
  normalise = numeric()
)

read_optimisation <- function(path, id, names = colnames(optimisation.table)){
  optimisation <- read.table(file = paste(path, "optimisation.csv", sep = "/"),
                             sep = ",", header = FALSE)
  dt <- data.table(matrix(c(id,as.numeric(optimisation[1,])), nrow = 1))
  colnames(dt) <- names
  return(dt)
}


plot_results <- function(data = data.exp.grouped,
                         path.analysis,
                         path.output = path.analysis,
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
        data.model.list = data.model.best.list[c("single", "receptors", as.character(id))],
        plot.title = paste(plot.title, "place:", i, ", id:", id, ",", sep = ""),
        filename = paste(path.analysis, plot.title, "_", i, ".pdf", sep = ""),
        plot.save = FALSE,
        plot.return = TRUE
      )
      
    })
  }
  pdf(file = paste(path.output, plot.title, ".pdf", sep = ""),
      width = 20,
      height = 12,
      useDingbats = FALSE)
  print(marrangeGrob(unlist(plot.list, recursive = FALSE), ncol = grid.ncol, nrow =  grid.nrow))
  dev.off()
}
#### ####
path.optimisation <- paste(path.output, "cmaes/mvn/2017-02-03/", sep = "/")
path.optimisation.data <- paste(path.optimisation, "data/", sep = "/")

stimulation.list.all <- ((data.exp %>%  distinct(stimulation))[-1])$stimulation

data.exp.grouped <-  read.table(
  file = paste(path.optimisation, "data_exp_grouped.csv", sep = ""),
  sep = ",",
  header = TRUE)

data.exp.grouped <- data.exp.grouped %>% group_by(priming, stimulation, time) %>% mutate(intensity_sd = var(intensity))

data.model.list <- list()
path.single <- paste(path.optimisation.data, "single", sep = "/")
optimisation.table <- rbind(optimisation.table,
                            read_optimisation(path = path.single,
                                              id = "single",
                                              names = colnames(optimisation.table)))
data.model.list[["single"]] <- read.table(
  file = paste(path.single, "data_model.csv", sep = "/"),
  sep = ",",
  header = TRUE)

path.receptors <- paste(path.optimisation.data, "receptors", sep = "/")
optimisation.table <- rbind(optimisation.table, read_optimisation(path = path.receptors, id = "receptors"))
data.model.list[["receptors"]] <- read.table(
  file = paste(path.receptors, "data_model.csv", sep = "/"),
  sep = ",",
  header = TRUE)
