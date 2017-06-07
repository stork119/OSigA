### ### ### ### ###
### jakstat_data 
### ### ### ### ###

#### libraries ####

#### input ####
path.resources <- "resources/"
path.data         <- paste(path.resources, "input/data/", sep = "")
path.output       <- paste(path.resources, "output/", sep = "")
dir.create(path.output, recursive = TRUE, showWarnings = TRUE)

path.data.output <- "resources/output/data/"
path.data.input  <- "resources/input/data"

#### ####

# data.exp <- data.table(position = character(),
#                        cell  = numeric(),
#                        priming = numeric(),
#                        stimulation = numeric(),
#                        time = numeric(),
#                        intensity = numeric())

#### KA10 ####
# data.exp.c <- read.table(file = paste(path.data, "KA10_C.csv", sep = "/"),
#                          header = TRUE,
#                          sep = "\t")
# 
# data.exp <- rbind(data.exp,
#                   data.table(position = data.exp.c$PositionName,
#                              cell = data.exp.c$ObjectNumber,
#                              priming = data.exp.c$group.1.1,
#                              stimulation = data.exp.c$group.2.1,
#                              time = data.exp.c$compare.1.1,
#                              intensity = data.exp.c$ShrinkedNuclei.Intensity))
# 
# data.exp.ifnb <- read.table(file = paste(path.data, "KA10_IFNB.csv", sep = "/"),
#                             header = TRUE,
#                             sep = "\t")
# data.exp <- rbind(data.exp,
#                   data.frame(position = data.exp.ifnb$PositionName,
#                              cell = data.exp.ifnb$ObjectNumber,
#                              priming = data.exp.ifnb$group.1.1,
#                              stimulation = data.exp.ifnb$group.2.1,
#                              time = data.exp.ifnb$compare.1.1,
#                              intensity = data.exp.ifnb$ShrinkedNuclei.Intensity))
# 
# data.exp$logintensity <- log(data.exp$intensity)
# data.exp.unique <- distinct(data.exp, priming, stimulation, time)
# data.exp.grouped <- data.exp %>% group_by(priming, stimulation, time)
# data.exp.grouped.all <- data.exp %>% group_by(priming, stimulation, time)
#### ####
data.list <- read_data(path = path.data.input)

data.list$data.exp %>% dplyr::distinct(file)

data.list$data.exp.norm <- get_equal_data(data = data.list$data.exp,
                                          sample_size = 1000)

data.list$data.exp.summarise <- 
  data.list$data.exp.norm %>% 
  dplyr::group_by(priming,
                  stimulation,
                  time) %>%
  summarise(m.norm = mean(intensity),
            sd.norm   = var(intensity))

data.list$data.exp.summarise<-
  data.list$data.exp.summarise %>%
  dplyr::mutate(mean.lmvn = lmvn.mean(m = m.norm, sd = sd.norm),
                sd.lmvn = lmvn.sd(m = m.norm, sd = sd.norm))

#### ####
#data.exp$time <- data.exp$time + 5
#tmesh.exp <- unique(data.exp$time)
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
# 
# data.exp.grouped.all$position
# 
# path <- paste(path.output, "data", sep = "/")
# dir.create(path)
# g <- ggplot(data.exp.grouped.all %>% filter(priming == 0), aes(x = factor(time - 5), y = intensity)) + 
#   facet_grid(~stimulation) + 
#   geom_boxplot() + 
#   theme_jetka() + ggtitle("MEFs IFNG stimulation") + ylim(c(0,1000))
# ggsave(filename = paste(path.output, "data", "control.pdf", sep = "/"), width = 24, height = 12, useDingbats = FALSE, plot = g  )
# 
# g <- ggplot(data.exp.grouped.all %>% filter(priming == 1000), aes(x = factor(time - 5), y = intensity)) + 
#   facet_grid(~stimulation) + 
#   geom_boxplot() + 
#   theme_jetka() + ggtitle("MEFs IFNB priming IFNG stimulation") + ylim(c(0,1000))
# ggsave(filename = paste(path.output, "data", "priming.pdf", sep = "/"), width = 24, height = 12, useDingbats = FALSE, plot = g  )
