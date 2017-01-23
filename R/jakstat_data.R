### ### ### ### ###
### jakstat_data 
### ### ### ### ###

path.resources <- "resources/"
path.data         <- paste(path.resources, "input/data/", sep = "")
path.output       <- paste(path.resources, "output/", sep = "")
dir.create(path.output, recursive = TRUE, showWarnings = TRUE)

#### ####

data.exp <- data.frame(type = character(), 
                       postion = character(),
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
                  data.frame(type = "C",
                             position = data.exp.c$PositionName,
                             cell = data.exp.c$ObjectNumber,
                             priming = data.exp.c$group.1.1,
                             stimulation = data.exp.c$group.2.1,
                             time = data.exp.c$compare.1.1,
                             intensity = data.exp.c$ShrinkedNuclei.Intensity))

data.exp.ifnb <- read.table(file = paste(path.data, "KA10_IFNB.csv", sep = "/"),
                            header = TRUE,
                            sep = "\t")
data.exp <- rbind(data.exp,
                  data.frame(type = "C",
                             position = data.exp.ifnb$PositionName,
                             cell = data.exp.ifnb$ObjectNumber,
                             priming = data.exp.ifnb$group.1.1,
                             stimulation = data.exp.ifnb$group.2.1,
                             time = data.exp.ifnb$compare.1.1,
                             intensity = data.exp.ifnb$ShrinkedNuclei.Intensity))


#### ####
data.exp$time <- data.exp$time + 5
tmesh.exp <- unique(data.exp$time)
data.exp$logintensity <- log(data.exp$intensity)

data.exp.background <- mean(data.exp[data.exp$time == 5,]$intensity)
