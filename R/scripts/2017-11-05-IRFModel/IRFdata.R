### ###
### IRF data
### ###

#### organize data IRF ####
irf_exp <- "2017-07-25-KS27" # "2017-07-18-KS25"
irfmodel.data.list <- list()

irfmodel.data.list$irf <- 
  poster.data.list[[irf_exp]] %>%
  dplyr::filter(time.2.1 == 120)

irfmodel.data.list$irf <- 
  do.call(
    rbind,
    list(
      irfmodel.data.list$irf,
      poster.data.list[[irf_exp]] %>%
        dplyr::filter(time.2.1 == 0) %>%
        dplyr::mutate(stimulation.1.1 = 0)
    )
  ) %>%
  dplyr::mutate(priming = priming.1.1) %>%
  dplyr::mutate(time = time.2.1) %>%
  dplyr::mutate(stimulation = stimulation.1.1) %>%
  dplyr::mutate(response = Intensity_MeanIntensity_Alexa) %>%
  dplyr::mutate(logresponse = log(response))

irfmodel.data.list$irf <- get_equal_data(data = irfmodel.data.list$irf)

irfmodel.data.list$irfsum <- 
  irfmodel.data.list$irf %>% 
  dplyr::group_by(stimulation) %>% 
  dplyr::summarise(
    logresponse.mean = mean(logresponse),
    logresponse.sd = sd(logresponse)
  ) %>%
  dplyr::mutate(
    logresponse = logresponse.mean
  )

#### organise data pSTAT ####
pstat_exp <-"2016-01-28-KA11-C" #"2016-01-26-KA10-C" ## "2016-01-28-KA11-C"
irfmodel.data.list$pSTAT <- 
  poster.data.list[[pstat_exp]] %>%
  dplyr::filter(time.1.1 == 30) 

irfmodel.data.list$pSTAT <- 
  do.call(
    rbind,
    list(
      irfmodel.data.list$pSTAT,
      poster.data.list[[pstat_exp]] %>%
        dplyr::filter(time.1.1 == 0) %>%
        dplyr::mutate(stimulation.1.1 = 0)
    )
  ) %>%
  dplyr::mutate(priming = priming.1.1) %>%
  dplyr::mutate(time = time.1.1) %>%
  dplyr::mutate(stimulation = stimulation.1.1) %>%
  dplyr::mutate(response = Intensity_MeanIntensity_Alexa) %>%
  dplyr::mutate(logresponse = log(response))

irfmodel.data.list$pSTAT <- get_equal_data(data = irfmodel.data.list$pSTAT)

irfmodel.data.list$pSTATsum <- 
  irfmodel.data.list$pSTAT %>% 
  dplyr::group_by(stimulation) %>% 
  dplyr::summarise(
    logresponse.mean = mean(logresponse),
    logresponse.sd = sd(logresponse)
  ) %>%
  dplyr::mutate(
    logresponse = logresponse.mean
  )
#### organize data all ####
data.raw.sum <-
  irfmodel.data.list$pSTATsum %>% 
  dplyr::full_join(irfmodel.data.list$irfsum, by = c("stimulation")) %>%
  dplyr::mutate(pstat = logresponse.x, 
                pstat.sd = logresponse.sd.x, 
                irf   = logresponse.y,
                irf.sd = logresponse.sd.y, 
                type  = "data") %>%
  dplyr::select(stimulation, pstat, pstat.sd, irf, irf.sd, type) %>%
  data.frame()

#### stimulation ####
stimulations.pSTAT <- irfmodel.data.list$pSTAT %>% dplyr::distinct(stimulation) 
stimulations.irf <- irfmodel.data.list$irf %>% dplyr::distinct(stimulation)
stimulations  <- sort((stimulations.pSTAT %>% dplyr::inner_join(stimulations.irf))$stimulation)
