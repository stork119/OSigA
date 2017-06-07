### ###
### model sensitivity analysis ###
### ###

#### ####
### ###
### 2017-06-01
### ###
source("~/Documents/ClustSense/ClusteringIdentifiability/R/SACluster.R")

delta            <- 0.95
zeta             <- 1
delimeter <- ","
directory.models <- '/media/knt/D/KN/ClustSense/Modelling/lib/StochSense/StochSens_Ver2.3/models/JAKSTAT_means/'
#### NFKB EXAMPLE ####
directory.folder <- paste(directory.models,
                          'Output/2017-05-30-priming/',
                          sep = "")

directory.SM     <- paste(directory.folder, 'sense.csv', sep = "")
#directory.FIM     <- paste(directory.folder, 'FIM.csv', sep = "")
directory.output <- paste(directory.folder,
                          "clustering/",
                          'delta/',
                          '/', 
                          sep = "")
directory.output.sens <- paste(directory.folder,
                          "sensitivity/",
                          '/', 
                          sep = "")

directory.params <- paste(directory.models,
                          'Input/parameters-names.csv',
                          sep = "")

dir.create(path = directory.output, recursive = TRUE, showWarnings = FALSE)

SM <- as.matrix(read.csv(directory.SM, sep = delimeter, header = FALSE))

params <- read.csv(directory.params, header = FALSE, sep = delimeter) ## HORIZONTAL

if( ncol(SM) != nrow(params)){
  print("Error. Wrong number of parameters")
}

x <- SMC(S = SM, 
         labels = t(params),
         names = t(params),
         zeta = zeta,
         delta = delta,
         dir.output =directory.output )

x$FIM_all
cluster <- clusterident.SMC(x)

directory.dendrogram <- paste(directory.output, "_dendrogram", ".pdf", sep = "")
pdf(directory.dendrogram, width = 12, height = 6)
plotDendogram.SMC(x,
                  cluster = cluster,
                  fig      = c(0,1,0,1),
                  new      = FALSE,
                  labels.args = list(cex = 0.5)
)
dev.off()

directory.fi <- paste(directory.output, "_fi", ".pdf", sep = "")
pdf(directory.fi, width = 6, height = 6)
p <- plotFI.SMC(x,
                cluster = cluster,
                fig      = c(0,1,0,1),
                new      = FALSE
)
dev.off()

directory.barplot <- paste(directory.output, "_barplot", ".pdf", sep = "")
pdf(directory.barplot, width = 12, height = 6)
barplot.SMC(x, cluster = cluster, fig = c(0,1,0,1))
dev.off()

directory.barplot <- paste(directory.output, "_barplot", ".pdf", sep = "")
pdf(directory.barplot, width = 12, height = 6)
barplot.SMC(x, cluster = cluster, fig = c(0,1,0,1))
dev.off()

directory.cov <- paste(directory.output, "_cov", ".pdf", sep = "")
pdf(directory.cov, width = 6, height = 6)
plotEV.Cov.SMC(x, cluster = cluster, fig = c(0,1,0,1))
dev.off()

directory.cor <- paste(directory.output, "_cor", ".pdf", sep = "")
pdf(directory.cor, width = 6, height = 6)
plotEV.Cor.SMC(x, cluster = cluster, fig = c(0,1,0,1))
dev.off()

directory.sa <- paste(directory.output, "_sa", ".csv", sep = "")
data.sa <- sa.SMC(x, cluster = cluster, header=FALSE)
write.csv(x = data.sa, file = directory.sa)

#### Sensitivity ####
df <- data.frame(params = params[,1], sens =  diag(x$FIM_all))
g.list[["sensitivity"]] <- 
  ggplot(data = df, mapping = aes(x = params, y = sens)) + geom_bar(stat = "identity") + 
    do.call(what = theme_jetka, args = plot.args) +
    ggtitle("Sensitivty analysis")
ggsave(filename = "sensitivity.pdf",
       path = directory.output.sens, 
       plot = g.list[["sensitivity"]],width  = plot.args$width, 
       height = plot.args$height, 
       useDingbats = plot.args$useDingbats)
#### observable STAT 1 only ####
###  Goal: which parameters, the concentration of STAT1 is the most sensitive to.

SM.stat1 <- SM[seq(from = 1, to = nrow(SM), by = 17),]
time.list <- seq(0,120,5)
time.list.ind <- which(time.list %in% tmesh[tmesh.list])
SM.stat1 <- SM.stat1[time.list.ind,]
SM.stat1 <- rbind(SM.stat1,SM.stat1)
SM.stat1 <- SM.stat1[,-5]
x.stat1 <- SMC(S = SM.stat1, 
         labels = t(params$V1[-5]),
         names = t(params$V1[-5]),
         zeta = zeta,
         delta = delta,
         dir.output =directory.output.sens )

# cluster <- clusterident.SMC(x.stat1)
# directory.dendrogram <- paste(directory.output.sens, "dendrogram_stat1", ".pdf", sep = "")
# pdf(directory.dendrogram, width = 12, height = 6)
# plotDendogram.SMC(x.stat1,
#                   cluster = cluster,
#                   fig      = c(0,1,0,1),
#                   new      = FALSE,
#                   labels.args = list(cex = 0.5)
# )
# dev.off()

df <- data.frame(params = params[-5,1], sens =  diag(x.stat1$FIM_all))
g.list[["sensitivity.stat1"]] <- 
  ggplot(data = df, mapping = aes(x = params, y = sens)) + geom_bar(stat = "identity") + 
  do.call(what = theme_jetka, args = plot.args) +
  ggtitle("Sensitivty analysis STAT1")
ggsave(filename = "sensitivity.stat1.pdf",
       path = directory.output.sens, 
       plot = g.list[["sensitivity.stat1"]],width  = plot.args$width, 
       height = plot.args$height, 
       useDingbats = plot.args$useDingbats)


#### primig vs single ####

return_SM_pstat1 <- function(path = '2017-05-30-priming', output.path = "pstat1/"){
directory.folder <- paste(directory.models,
                          'Output',
                          path,
                          "/",
                          sep = "/")

directory.SM     <- paste(directory.folder, 'sense.csv', sep = "")
#directory.FIM     <- paste(directory.folder, 'FIM.csv', sep = "")
directory.output <- paste(directory.folder,
                          output.path,
                          '/', 
                          sep = "")


directory.params <- paste(directory.models,
                          'Input/parameters-names.csv',
                          sep = "")
directory.dendrogram <- paste(directory.output, "dendrogram", ".pdf", sep = "")
dir.create(path = directory.output, recursive = TRUE, showWarnings = FALSE)

params <- read.csv(directory.params, header = FALSE, sep = delimeter)[-5,] ## HORIZONTAL

SM.priming <- as.matrix(read.csv(directory.SM, sep = delimeter, header = FALSE))

SM.priming <- SM.priming[,-5]
SM.priming.pstat1 <- SM.priming[seq(from = 13, to = nrow(SM.priming), by = 17),]
time.list <- seq(0,120,5)
time.list.ind <- which(time.list %in% tmesh[tmesh.list])
SM.priming.pstat1 <- SM.priming.pstat1[time.list.ind,]
SM.priming.pstat1 <- do.call(rbind, list(SM.priming.pstat1,SM.priming.pstat1,SM.priming.pstat1,SM.priming.pstat1,SM.priming.pstat1))
SM.priming.pstat1 <- do.call(rbind, list(SM.priming.pstat1,SM.priming.pstat1,SM.priming.pstat1,SM.priming.pstat1,SM.priming.pstat1))
SM.priming.pstat1 <- do.call(rbind, list(SM.priming.pstat1,SM.priming.pstat1,SM.priming.pstat1,SM.priming.pstat1,SM.priming.pstat1))


x.priming.pstat1 <- SMC(S = SM.priming.pstat1, 
               labels = params,
               names = params,
               zeta = zeta,
               delta = delta,
               dir.output = directory.output )

df <- data.frame(params = params, sens =  diag(x.priming.pstat1$FIM_all))

g.list<- list()
g.list[["sensitivity"]] <- 
  ggplot(data = df, mapping = aes(x = params, y = sens)) + geom_bar(stat = "identity") + 
  do.call(what = theme_jetka, args = plot.args) +
  ggtitle("Sensitivty analysis pSTAT1")
ggsave(filename = "sensitivity-pstat.pdf",
       path = directory.output, 
       plot = g.list[["sensitivity"]],width  = plot.args$width, 
       height = plot.args$height, 
       useDingbats = plot.args$useDingbats)

cluster.priming.pstat1 <- clusterident.SMC(x.priming.pstat1)


pdf(directory.dendrogram, width = 12, height = 6)
print(plotDendogram.SMC(x.priming.pstat1,
                  cluster = cluster.priming.pstat1,
                  fig      = c(0,1,0,1),
                  new      = FALSE,
                  labels.args = list(cex = 0.5)
))
dev.off()
return(SM.priming.pstat1)
}
##### ####
SM.list <-list(priming = return_SM_pstat1(path = '2017-05-30-priming'),
               single = return_SM_pstat1(path = '2017-05-30-single'))


SM <- rbind(SM.list$priming,SM.list$single)
path = '2017-05-30'
directory.folder <- paste(directory.models,
                          'Output',
                          path,
                          "/",
                          sep = "/")
directory.output <- paste(directory.folder,
                          output.path,
                          '/', 
                          sep = "")
directory.dendrogram <- paste(directory.output, "dendrogram", ".pdf", sep = "")

dir.create(directory.output, recursive = TRUE)

x.pstat1 <- SMC(S = SM, 
                        labels = params,
                        names = params,
                        zeta = zeta,
                        delta = delta,
                        dir.output = directory.output )

df <- data.frame(params = params, sens =  diag(x.pstat1$FIM_all))

g.list<- list()
g.list[["sensitivity"]] <- 
  ggplot(data = df, mapping = aes(x = params, y = sens)) + geom_bar(stat = "identity") + 
  do.call(what = theme_jetka, args = plot.args) +
  ggtitle("Sensitivty analysis pSTAT1")
ggsave(filename = "sensitivity-pstat.pdf",
       path = directory.output, 
       plot = g.list[["sensitivity"]],width  = plot.args$width, 
       height = plot.args$height, 
       useDingbats = plot.args$useDingbats)

cluster.pstat1 <- clusterident.SMC(x.pstat1)


pdf(directory.dendrogram, width = 12, height = 6)
print(plotDendogram.SMC(x.pstat1,
                        cluster = cluster.pstat1,
                        fig      = c(0,1,0,1),
                        new      = FALSE,
                        labels.args = list(cex = 0.5)
))
dev.off()
