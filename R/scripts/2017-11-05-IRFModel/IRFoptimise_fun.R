### ###
### IRF optimise_fun
### ###

#### optimisation function ####
optimise.fun <- function(par,
                             stimulations,
                             data.raw.list,
                             ranges.factor,
                             ranges.base,
                             ranges.opt,
                             model_fun,
                             data.model.colnames,
                             ...
){
  tryCatch({
  params <- ranges.factor
  params[ranges.opt] <- ranges.factor[ranges.opt]*ranges.base[ranges.opt]^par
  
  data.model <- model_fun(
    stimulations = stimulations,
    params = params,
    par = par,
    #data.model.ps1 = data.model.ps1)
    #ranges = ranges.irf1, data.sample.ps1 = data.sample.sdnonconst.ps1)
    #
    ...
  )
  likelihood.list <- foreach( data.i = 1:length(data.raw.list) ) %do% {
    data  <- data.raw.list[[data.i]]
    normalise <- (data %>%
                    dplyr::mutate(normalise = (logresponse^2)/(logresponse.sd^2)) %>%
                    dplyr::summarise(normalise = mean(normalise)))$normalise
    likelihood <-
      (data %>% 
         left_join(
           by =  "stimulation",
           ((data.model %>% 
               dplyr::mutate_(
                 "logmodel" = data.model.colnames[data.i],
                 "logmodel.sd" = paste(data.model.colnames[data.i], "sd", sep = ".")
               ))[, 
                  c("logmodel", 
                    "logmodel.sd",
                    "stimulation")])) %>%
         #dplyr::mutate(likelihood = ((logresponse - logmodel)^2)/(2*(logmodel.sd^2))  + log(logmodel.sd)) %>%
         dplyr::mutate(likelihood = ((logresponse - logmodel)^2)/(logresponse.sd^2)) %>%
         #  dplyr::mutate(likelihood = ((logresponse - logmodel)^2)/(logresponse^2)) %>%
         dplyr::summarise(likelihood = sum(likelihood)))$likelihood
    return(likelihood/normalise)         
  }
  likelihood <- do.call(what = sum, args = likelihood.list)
  # print(likelihood)
  # print("\n")
  # print(par)
  return(as.numeric(likelihood))
  }, error = function(e){print(e)})
  return(Inf)
}