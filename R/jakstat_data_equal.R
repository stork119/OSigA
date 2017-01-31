###
### jakstat_data_equal ####
###


get_equal_data  <- function(data){

d.distinct <- data %>% distinct(priming, stimulation, time) 
d.i <- d.distinct[1,]
data.equal <- data %>% 
  filter(stimulation == d.i$stimulation, priming == d.i$priming, time == d.i$time) %>% 
  sample_n(1000, replace = TRUE)

for(i in 2:nrow(d.distinct)){
  d.i <- d.distinct[i,]
  data.equal <- data.equal %>%
    rbind(data %>% 
            filter(stimulation == d.i$stimulation, priming == d.i$priming, time == d.i$time) %>% 
            sample_n(1000, replace = TRUE))
}
 return(data.equal)
}