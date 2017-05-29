### ### ### 
### jakstat_data_equal ###
### ### ###

#### ####
### function returns data with equal number of elements in all subrgroups (priming, stimulation, time)
get_equal_data  <- function(
  data,
  sample_size = 1000){
  
  d.distinct <- data %>%
    distinct(priming,
             stimulation,
             time) 
  d.i <- d.distinct[1,]
  data.equal <- data %>% 
    filter(stimulation == d.i$stimulation,
           priming == d.i$priming,
           time == d.i$time) %>% 
    sample_n(sample_size, replace = TRUE)
  
  for(i in 2:nrow(d.distinct)){
    d.i <- d.distinct[i,]
    data.equal <- data.equal %>%
      rbind(data %>% 
              filter(stimulation == d.i$stimulation,
                     priming == d.i$priming,
                     time == d.i$time) %>% 
              sample_n(sample_size, replace = TRUE))
  }
   return(data.equal)
}