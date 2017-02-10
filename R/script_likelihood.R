setwd("~/Documents/modelling/")
source("R/parallel_computing.R")

path.optimisation <- paste(path.output, "cmaes/normalized/2017-02-04-5/", sep = "/")
path.optimisation.data <- paste(path.optimisation, "data/", sep = "/")

ids.list <- list.files(path.optimisation.data)
ids.list <- optimisation.table$id

st.likelihood <- data.table(id= character(),
                            priming = numeric(),
                            control = numeric(),
                            all = numeric())
g.list <- list()
ids.list.which <- which(ids.list %in% st.likelihood$id)
if(length(ids.list.which) != 0){
  ids.list <- ids.list[-ids.list.which]
}


for(id in ids.list){

filename.data_model <- "data_model.csv"
data.exp.grouped <-  read.table(
  file = paste(path.optimisation, "data_exp_grouped.csv", sep = ""),
  sep = ",",
  header = TRUE)

data.model <- read.table(
  file = paste(path.optimisation.data, id, filename.data_model, sep = "/"),
  sep = ",",
  header = TRUE)

result <-
  likelihood(
    data.model = data.model, 
    data.exp.grouped = data.exp.grouped,
    fun.likelihood = fun.normalised)

data.model$likelihood <- result

st.likelihood <- rbind(st.likelihood,
                       data.table(id= id,
                                  priming = as.numeric(data.model %>% filter(priming == 1000) %>% summarise(count = sum(likelihood))),
                                  control = as.numeric(data.model %>% filter(priming == 0) %>% summarise(count = sum(likelihood))),
                                  all = as.numeric(data.model %>% summarise(count = sum(likelihood))) )
                       )

g.list[[id]] <- ggplot(data.model, aes(x = factor(time - 5), y = likelihood)) +
  geom_point() +
  facet_grid(priming~stimulation) + 
  theme_jetka() +
  ylim(c(0,0.3)) +
  ggtitle(id)
}

st.likelihood <- st.likelihood[order(st.likelihood$all),]

# pdf(file = paste(path.optimisation, "likelihood", sep = "/"),
#     widt= 24,
#     height =12,
#     useDingbats = FALSE)
# print(g.list)
# dev.off()


pdf(file = paste(path.optimisation, "likelihood_priming", sep = "/"),
    widt= 24,
    height =12,
    useDingbats = FALSE)
for(id in st.likelihood[order(st.likelihood$priming),]$id[1:25]){
print(g.list[[as.character(id)]])
}
dev.off()