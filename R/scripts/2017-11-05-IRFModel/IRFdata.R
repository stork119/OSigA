### ###
### IRF data
### ###

#### organize data IRF ####
irfmodel.data.list <- list()

irfmodel.data.list$irf <- 
  poster.data.list$`Gamma-IRF-Nuclei` %>%
  dplyr::filter(time == 120)
irfmodel.data.list$irf <- 
  do.call(
    rbind,
    list(
      irfmodel.data.list$irf,
      poster.data.list$`Gamma-IRF-Nuclei` %>%
        dplyr::filter(time == 0) %>%
        dplyr::mutate(stimulation = 0)
    )
  ) %>%
  dplyr::mutate(response = Intensity_MeanIntensity_Alexa) %>%
  dplyr::mutate(logresponse = log(response))
irfmodel.data.list$irf <- get_equal_data(irfmodel.data.list$irf)

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
irfmodel.data.list$pSTAT <- 
  poster.data.list$`Gamma-pStat-Nuclei` %>%
  dplyr::filter(time == 30) %>%
  dplyr::mutate(response = Intensity_MeanIntensity_Alexa) %>%
  dplyr::mutate(logresponse = log(response))

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
irfmodel.data.list$pSTAT <- get_equal_data(irfmodel.data.list$pSTAT)
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

