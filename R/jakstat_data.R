### ### ### ### ###
### jakstat_data 
### ### ### ### ###

#### libraries ####

#### input ####
path.resources <- "resources/"
path.data         <- paste(path.resources, "input/data/", sep = "")
path.output       <- paste(path.resources, "output/", sep = "")
dir.create(path.output, recursive = TRUE, showWarnings = TRUE)

#### ####

data.exp <- data.table(position = character(),
                       cell  = numeric(),
                       priming = numeric(),
                       stimulation = numeric(),
                       time = numeric(),
                       intensity = numeric())

#### KA10 ####
data.exp.c <- read.table(file = paste(path.data, "KA10_C.csv", sep = "/"),
                         header = TRUE,
                         sep = "\t")

data.exp <- rbind(data.exp,
                  data.table(position = data.exp.c$PositionName,
                             cell = data.exp.c$ObjectNumber,
                             priming = data.exp.c$group.1.1,
                             stimulation = data.exp.c$group.2.1,
                             time = data.exp.c$compare.1.1 + 5,
                             intensity = data.exp.c$ShrinkedNuclei.Intensity))

data.exp.ifnb <- read.table(file = paste(path.data, "KA10_IFNB.csv", sep = "/"),
                            header = TRUE,
                            sep = "\t")
data.exp <- rbind(data.exp,
                  data.frame(position = data.exp.ifnb$PositionName,
                             cell = data.exp.ifnb$ObjectNumber,
                             priming = data.exp.ifnb$group.1.1,
                             stimulation = data.exp.ifnb$group.2.1,
                             time = data.exp.ifnb$compare.1.1 + 5,
                             intensity = data.exp.ifnb$ShrinkedNuclei.Intensity))


#### ####
#data.exp$time <- data.exp$time + 5
# tmesh.exp <- unique(data.exp$time)
# data.exp$logintensity <- log(data.exp$intensity)
# 
# data.exp.background <- mean(data.exp[data.exp$time == 5,]$intensity)
# 
# data.exp.unique <- distinct(data.exp, priming, stimulation, time)
# 
# 
# data.exp.grouped <- data.exp %>% group_by(priming, stimulation, time)
# data.mean <- dplyr::summarise(data.exp.grouped, mean = mean(intensity), logmean = mean(logintensity))
# dplyr::filter(data.exp.grouped, priming == 0, time == 5, stimulation ==0.01)

#ggplot(data.mean, mapping = aes(x = time, y = logmean, color = factor(priming), group = interaction(priming, stimulation))) + geom_line()
