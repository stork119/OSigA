###
### jakstat_data_equal ####
###

d <- data.exp.grouped %>% filter(stimulation %in% stimulation.list)
d.distinct <- d %>% distinct(priming, stimulation, time) 

d.i <- d.distinct[1,]
data.exp.grouped.equal <- data.exp.grouped %>% 
  filter(stimulation == d.i$stimulation, priming == d.i$priming, time == d.i$time) %>% 
  sample_n(1000, replace = TRUE)

for(i in 2:nrow(d.distinct)){
  d.i <- d.distinct[i,]
  data.exp.grouped.equal <- data.exp.grouped.equal %>%
    rbind(data.exp.grouped %>% 
            filter(stimulation == d.i$stimulation, priming == d.i$priming, time == d.i$time) %>% 
            sample_n(1000, replace = TRUE))
}

write.table(x = data.exp.grouped.equal,
            file = paste(path.optimisation, "data_exp_grouped.csv", sep = ""),
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)