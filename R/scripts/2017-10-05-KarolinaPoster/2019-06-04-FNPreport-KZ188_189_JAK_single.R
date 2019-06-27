### ### ###
### JAK  KZ 188 189
### ### 
library(dplyr)
library(SysBioSigTheme)
library(ITRC)
#### computeInhibitorsFractions ####
computeInhibitorsFractions <- 
  function(
    data.curve,
    data.curve.distinct,
    response.type,
    levels.list,
    level.min,
    level.max,
    bootstrap = TRUE,
    bootstrap.number = 8,
    bootstrap.sample_size = 1000,
    bootstrap.test.sample = TRUE,
    bootstrap.test.number = 4,
    parallel_cores = 8,
    ...
  ){
    
    foreach(level_ = levels.list) %do% {
      
      itrc.model <- ITRC::ITRC(
        data = data.curve %>%
          dplyr::filter(level %in% c(level.min, level.max, level_)), 
        signal = "level",
        response = response.type,
        bootstrap.number = bootstrap.number,
        bootstrap = bootstrap,
        bootstrap.sample_size = bootstrap.sample_size,
        bootstrap.test.sample = bootstrap.test.sample,
        bootstrap.test.number = bootstrap.test.number,
        parallel_cores = parallel_cores
      )
      itrc.model$confusion.table %>% 
        dplyr::filter(max.signal == level.max) %>%
        dplyr::filter(level == level_) %>%
        dplyr::left_join(
          (data.frame(
            class = c(level.min, level_, level.max),
            type = c("complete", "partial", "no_inhibition")
          )),
          by = "class"
        ) %>% dplyr::ungroup() ->
        confusion.table
      confusion.table %>% 
        dplyr::filter(type != "partial") %>% 
        rbind(.,
              confusion.table %>% 
                dplyr::filter(type == "no_inhibition") %>%
                dplyr::mutate(prob =  1 - prob) %>% 
                dplyr::mutate(type = "partial")
        )
      
    } %>% 
      do.call(what = rbind,
              args = .) %>%
      dplyr::left_join(data.curve.distinct,
                       by = "level") ->
      inhibitor.df
  }
#### computeInhibitorsFractionsScalesFixed ####
computeInhibitorsFractionsScalesFixed <- 
  function(
    data.curve,
    data.curve.distinct,
    response.type,
    levels.list,
    level.min,
    level.max,
    bootstrap = TRUE,
    bootstrap.number = 8,
    bootstrap.sample_size = 1000,
    bootstrap.test.sample = TRUE,
    bootstrap.test.number = 4,
    parallel_cores = 8,
    quantile.prob = 0.1,
    cc_maxit = 100,
    lr_maxit = 1000,
    MaxNWts = 5000,
    ...
  ){
    
    df.res <-
      GetLogRegParameters(
        data = 
          data.curve %>%
          dplyr::filter(level %in% c(level.min, level.max)),
        response = response.type,
        signal = "level",
        signal.list = c(level.min, level.max))
    
    
    lr_model <- nnet::multinom(formula = df.res$formula_string,
                               data = df.res$data,
                               na.action = na.omit,
                               maxit = lr_maxit,
                               MaxNWts = MaxNWts)#,  model = FALSE )
    
    
    lr.fit <-
      predict(object  = lr_model,
              newdata = df.res$data)
    
    df.res$data$class <- as.numeric(as.character(lr.fit))
    
    df.res$data %>% 
      dplyr::rename(level = X) %>% 
      dplyr::rename_(.dots = setNames(
        object =df.res$response,
        nm = response.type)) %>%
      dplyr::filter(level == class) ->
      df.res$data
    
    fraction.min.bound <-
      quantile((df.res$data %>%
                  dplyr::filter(level %in% c(level.min)))[[response.type]],
               probs = 1 - quantile.prob)[[1]]
    
    fraction.max.bound <-
      quantile((df.res$data %>%
                  dplyr::filter(level %in% c(level.max)))[[response.type]],
               probs = quantile.prob)[[1]]
    
    bound.df <-
      data.frame(
        type = c("complete", "partial"),
        bound = c(fraction.min.bound,fraction.max.bound),
        stringsAsFactors = FALSE
      )
    foreach(level_ = levels.list) %do% {
      
      data.curve %>%
        dplyr::filter(level %in% c( level_)) %>%
        dplyr::mutate_(complete = paste(response.type, "<", fraction.min.bound)) %>%
        dplyr::mutate_(partial = paste(response.type, "<", fraction.max.bound)) %>%
        dplyr::group_by(level) %>%
        dplyr::summarise(partial = sum(partial)/n(),
                         complete = sum(complete)/n()) %>% 
        reshape2::melt(
          id.vars = c( "level"),
          variable.name = "type",
          value.name = "prob", factorsAsStrings = TRUE
        ) %>%
        dplyr::left_join(
          bound.df,
          by = "type"
        ) 
      
    } %>% 
      do.call(what = rbind,
              args = .) %>%
      dplyr::left_join(data.curve.distinct,
                       by = "level") ->
      inhibitor.df
    return(inhibitor.df)
  }
#### ####
#poster.path.list <- list()
#poster.path.list$inhibitor.input.dir <- paste(poster.path.list$input.date.dir, "JAK_inhibitors", sep = "/")
poster.path.list$inhibitor.input.dir <- paste("/media/knt/sdb2/KN/IFNsignalling/2019-06-17-JAKsignle/resources/", sep = "/")



folder1 <- ""

folders.to.remove <- c(paste0("G", c(paste0("0", 1:9),10:12)), paste0("H", c(paste0("0", 1:9),10:12)))

## KZ188
data.KZ188 <- (read.table(paste(poster.path.list$inhibitor.input.dir, folder1, "ShrinkedNuclei_KZ188_pS1.csv", sep="/"), 
                          header=TRUE, sep=",") %>% 
                 SysBioSigTheme::normalize_data(data = .))$data
data.KZ188 %>% dplyr::filter(!(well.name %in% folders.to.remove)) -> data.KZ188

data <- data.KZ188
title_ <- "KZ188"
response.df <-
  data.frame( 
    type = c("Intensity_MeanIntensity_Alexa"),
    label = c("pS1"),stringsAsFactors = FALSE)

## KZ189
data.KZ189 <- (read.table(paste(poster.path.list$inhibitor.input.dir, folder1, "ShrinkedNuclei_KZ189_pS3.csv", sep="/"), 
                          header=TRUE, sep=",") %>% 
                 SysBioSigTheme::normalize_data(data = .))$data

data.KZ189 %>% dplyr::filter(!(well.name %in% folders.to.remove)) -> data.KZ189
data.KZ189$stimulation.1.2 <- 0
data <- data.KZ189
title_ <- "KZ189"


response.df <-
  data.frame( 
    type = c("Intensity_MeanIntensity_Alexa"),
    label = c("pS3"),stringsAsFactors = FALSE)

## KZ192
data.KZ191 <- (read.table(paste(poster.path.list$inhibitor.input.dir, folder1, "ShrinkedNuclei_KZ191_pS1pS3.csv", sep="/"), 
                          header=TRUE, sep=",") %>% 
                 SysBioSigTheme::normalize_data(data = .))$data
title_ <- "KZ191"
data <- data.KZ191

## KZ192
data.KZ192 <- (read.table(paste(poster.path.list$inhibitor.input.dir, folder1, "ShrinkedNuclei_KZ192_pS1pS3.csv", sep="/"), 
                          header=TRUE, sep=",") %>% 
                 SysBioSigTheme::normalize_data(data = .))$data
title_ <- "KZ192"
data <- data.KZ192
response.df <-
  data.frame( 
    type = c("Intensity_MeanIntensity_Alexa488",
             "Intensity_MeanIntensity_Alexa555"),
    label = c("pS1",
              "pS3"),stringsAsFactors = FALSE)





#### computation ####
inhibitor.list <- c(
  "inhibitor.1.1",
  "inhibitor.1.2",
  "inhibitor.1.3")

inhibitors.labels <- c("JAK1", "JAK2", "TYK2")


exp.df <-
  expand.grid(
    inhibitor.i = 1:length(inhibitors.labels),
    response.i = 1:nrow(response.df)
  )

foreach(exp.i = 1:nrow(exp.df)) %do% {
  response.i  <- exp.df[exp.i,]$response.i
  inhibitor.i <- exp.df[exp.i,]$inhibitor.i
  
  inhibitor_ = inhibitor.list[inhibitor.i]
  inhibitor.label <- inhibitors.labels[inhibitor.i]
  inhibitor.x <- inhibitor.list[inhibitor.list != inhibitor_][1]
  inhibitor.x.label <- inhibitors.labels[-inhibitor.i][1]
  inhibitor.y <- inhibitor.list[inhibitor.list != inhibitor_][2]
  inhibitor.y.label <- inhibitors.labels[-inhibitor.i][2]
  
  response.type <- response.df[response.i,"type"]
  response.label <- response.df[response.i,"label"]
  
  data %>% dplyr::filter_(
    paste(inhibitor.x, "== 0", "&", inhibitor.y, " == 0")) ->
    data.curve
  
  data.curve %>%
    dplyr::distinct_("stimulation.1.1", "stimulation.1.2", inhibitor_) ->
    data.curve.distinct
  
  data.curve.distinct %>%
    dplyr::filter_(paste(inhibitor_, "== 0")) %>% 
    dplyr::filter_("stimulation.1.1 == 0 & stimulation.1.2 == 0") %>%
    rbind(.,
          (data.curve.distinct %>%
             dplyr::filter_(paste(inhibitor_, "!= 0"))) %>% 
            dplyr::arrange_(inhibitor_)) %>%
    rbind(.,
          (data.curve.distinct %>%
             dplyr::filter_(paste(inhibitor_, "== 0")) %>% 
             dplyr::filter_("stimulation.1.1 != 0 | stimulation.1.2 != 0"))) ->
    data.curve.distinct
  data.curve.distinct$level <- 1:nrow(data.curve.distinct)
  
  data.curve %>%
    dplyr::left_join(data.curve.distinct,
                     by = c("stimulation.1.1", "stimulation.1.2", inhibitor_)) ->
    data.curve
  
  
  computeInhibitorsFractionsScalesFixed(
    data.curve = data.curve,
    data.curve.distinct = data.curve.distinct,
    response.type = response.type,
    levels.list   = (data.curve %>% 
                       dplyr::filter(level != min(level),
                                     level != max(level)
                       ) %>% 
                       dplyr::arrange(level) %>%
                       dplyr::distinct(level))[["level"]],
    level.min  = min(data.curve$level),
    level.max  = max(data.curve$level)) ->
    inhibitor.df
  inhibitor.df$inhibitor.i <- inhibitor.i
  inhibitor.df$inhibitor.label <- inhibitor.label
  inhibitor.df$response.type <- response.type
  inhibitor.df$response.label <-response.label
  
  return(inhibitor.df)
  
} ->
  inhibitor.df.list

foreach(inhibitor.df = inhibitor.df.list) %do% {
  colnames(inhibitor.df)[7] <- "inhibitor"
  return(inhibitor.df)
}  %>%
  do.call(what = rbind,args = .) ->
  inhibitor.df


ggplot(data = inhibitor.df %>% dplyr::filter(type != "no_inhibition"),
       aes_string(
         x = paste("factor(inhibitor)"),
         y = "prob",
         # ymin = "prob - prob.sd",
         # ymax = "prob + prob.sd",
         group = "type",
         color = "type")) +
  geom_line() +
  # geom_errorbar() +
  coord_cartesian(ylim = c(0,1)) +
  SysBioSigTheme::theme_sysbiosig() +
  facet_grid(inhibitor.label ~ response.label) +
  scale_color_viridis(discrete = TRUE, end = 0.8) ->
  g.scrc 

ggsave(filename = paste(poster.path.list$output.fnp.dir, 
                        paste0("inhibitors-distinctcurve-fixed-scrc-", title_,".pdf"), sep = "/"), 
       plot = g.scrc,
       width = 12,
       height = 12, useDingbats = FALSE
)
saveRDS(object = inhibitor.df,
        paste(poster.path.list$output.fnp.dir, 
              paste("inhibitors-scalefixed-scrc", title_, ".rds"),
              sep = "/"
        ))
#### ####

# inhibitor.df.res <- inhibitor.df
# ggplot(data.curve, aes_string(x = paste("log(",as.character(response.type),")"), color = "factor(level)", group = "level")) + geom_density() + SysBioSigTheme::theme_sysbiosig()
# 
# inhibitor.df %>% dplyr::filter(inhibitor.label == "TYK2", response.label0 == "pS3")

inhibitor.df %>% dplyr::filter(type != "no_inhibition") %>% 
  dplyr::select(level, type, inhibitor, inhibitor.label, response.label, prob) %>%
  dplyr::filter(response.label == "pS1") %>%
  dplyr::left_join(
    inhibitor.df %>% dplyr::filter(type != "no_inhibition") %>% 
      dplyr::select(level, type, inhibitor, inhibitor.label, response.label, prob) %>%
      dplyr::filter(response.label == "pS3"),
    by = c("level", "type", "inhibitor", "inhibitor.label" ), suffix = c(".pS1", ".pS3")
  ) %>%
  dplyr::mutate(prob.diff = prob.pS1 - prob.pS3) ->
  inhibitors.compared


ggplot(data = inhibitors.compared,
       aes_string(
         x = paste("factor(inhibitor)"),
         y = "prob.diff",
         group = "type")) +
  geom_line() +
  # geom_errorbar() +
  coord_cartesian(ylim = c(-1,1)) +
  SysBioSigTheme::theme_sysbiosig() +
  facet_grid(inhibitor.label ~ type) +
  scale_color_viridis(discrete = TRUE, end = 0.8) ->
  g.scrc.compare 
ggsave(filename = paste(poster.path.list$output.fnp.dir, 
                        paste0("inhibitors-distinctcurve-fixed-scrc-comparison",title_,".pdf"), sep = "/"), 
       plot = g.scrc.compare,
       width = 12,
       height = 12, useDingbats = FALSE
)

ggplot(data %>% 
  dplyr::filter(inhibitor.1.2 == 0, inhibitor.1.3 == 0) ,
  aes(x = log(Intensity_MeanIntensity_Alexa555), group = interaction(inhibitor.1.1, stimulation.1.1), color = interaction(inhibitor.1.1, stimulation.1.1))) +
    geom_density() +SysBioSigTheme::theme_sysbiosig()