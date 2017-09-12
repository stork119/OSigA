### ###
### likelihood
### ###

likelihood <- function(data.model,
                       data.exp.grouped,
                       data.exp.summarise,
                       fun.likelihood,
                       ...
){
  res <- sapply(1:nrow(data.model),
                function(data.model.i){
                  data.model.tmp <- data.model[data.model.i,]
                  zosia <-data.exp.grouped %>%
                    dplyr::filter(priming == data.model.tmp$priming,
                                  time == data.model.tmp$time,
                                  stimulation == data.model.tmp$stimulation) 
                  return((zosia%>%
                            dplyr::mutate(likelihood = 
                                            do.call(fun.likelihood,list(logintensity = logintensity, 
                                                                        intensity = intensity, 
                                                                        data.model.tmp = data.model.tmp,
                                                                        data.exp.summarise = data.exp.summarise))) %>%
                            summarise(likelihood.sum = 
                                        sum(likelihood)))$likelihood.sum)
                }
  )
  
  return(res)
}