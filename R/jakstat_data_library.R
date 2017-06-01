### ### ### 
### jakstat_data_equal ###
### ### ###

#### data normalization ####
normalization <- function(data.model,
                          m.scale = 400,
                          sd.scale = m.scale^2,
                          background,
                          epsilon = 1){
  data.model$m.norm <- data.model$m/m.scale + background
  data.model$sd.norm <- data.model$sd/sd.scale
  data.model$sd.norm <- sapply(data.model$sd.norm, function(sd.norm){ifelse(sd.norm < epsilon, epsilon, sd.norm)})
  data.model$time <- data.model$time - 5
  return(data.model)
}

#### get data equal ####
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
