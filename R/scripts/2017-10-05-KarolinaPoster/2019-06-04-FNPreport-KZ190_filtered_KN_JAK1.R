library(dplyr)
library(SysBioSigTheme)
library(ITRC)

poster.path.list$inhibitor.input.dir <- paste(poster.path.list$input.date.dir, "JAK_inhibitors", sep = "/")


folder1 <- ""
KZ190 <- read.table(paste(poster.path.list$inhibitor.input.dir, folder1, "ShrinkedNuclei.csv", sep="/"), 
                  header=TRUE, sep=",")

KZ190n <- SysBioSigTheme::normalize_data(KZ190)$data

KZ190n_filtered <- subset(KZ190n, ! KZ190n$well.name %in% c("A02", "B02", "B03",
        "B10", "B13", "C10", "C12", "C13", "D02", "D09", "D12", "E06", "E13",
        "F01", "F12", "A17", "A18", "B17", "B18", "C17", "C18", "D17", "D18"))

KZ190n_pS1 <- subset(KZ190n_filtered, KZ190n_filtered$antibody.1.1=="1:100")
KZ190n_pS3 <- subset(KZ190n_filtered, KZ190n_filtered$antibody.1.2=="1:200")

g.inhibitors <- list()
#### JAK2=0 ####

# color.limits <-
  # quantile(x = KZ190n_pS1[["Intensity_MeanIntensity_Alexa488"]],
        # probs = c(0.1,0.9))
# zakres: color.limits <- c(15, 300)
# mean-max: color.limits <- NULL
# color.limits <- c(15, max(KZ190n_pS1.sum$response.mean))

KZ190n_pS3_J2_0 <- subset(KZ190n_pS3, KZ190n_pS3$inhibitor.1.2=="0")

KZ190n_pS3_J2_0 %>% head()

KZ190n_pS3_J2_0 %>%
  dplyr::group_by(
    inhibitor.3.1,
    inhibitor.3.3) %>%
  dplyr::summarise(
    response.mean = mean(Intensity_MeanIntensity_Alexa555),
    response.median = median(Intensity_MeanIntensity_Alexa555)
  ) ->
  KZ190n_pS3_J2_0.sum


cbind(KZ190n_pS3_J2_0.sum %>% dplyr::ungroup(),
      KZ190n_pS3_J2_0.sum %>% dplyr::filter(
        inhibitor.3.1 == 0,
        inhibitor.3.3 == 0) %>% 
        dplyr::ungroup() %>%
        dplyr::rename(scale.mean = response.mean, 
                      scale.median = response.median) %>%
        dplyr::select(scale.mean, scale.median) %>% 
        dplyr::ungroup() ) %>%
  dplyr::mutate(
    relative.mean = response.mean/scale.mean,
    relative.median = response.median/scale.median
  ) ->
  KZ190n_pS3_J2_0.sum

color.limits <- c(min(KZ190n_pS3_J2_0.sum$response.median), max(KZ190n_pS3_J2_0.sum$response.median))
color.limits <- c(0, 1)
g.inhibitors[["JAK2-ps3"]] <-
  ggplot(data = KZ190n_pS3_J2_0.sum,
       mapping = aes(x = factor(inhibitor.3.1),
                     y = factor(inhibitor.3.3),
                     fill = relative.median,
                     label = round(relative.median,2)
       )) +
  geom_tile(color = "white") +
  viridis::scale_fill_viridis(
    name = "Relative response\npY-STAT3",
    limits = color.limits) + 
  ggtitle("JAK2 = 0")+
  xlab("JAK1")+
  ylab("TYK2")+
  ITRC::theme_itrc() + 
  geom_text(size = 8)
g.inhibitors[["JAK2-ps3"]]


#### TYK2=0 ####
KZ190n_pS1_T2_0 <- subset(KZ190n_pS1, KZ190n_pS1$inhibitor.1.3=="0")

KZ190n_pS1_T2_0 %>% head()

KZ190n_pS1_T2_0 %>%
  dplyr::group_by(
    inhibitor.3.1,
    inhibitor.3.2) %>%
  dplyr::summarise(
    response.mean = mean(Intensity_MeanIntensity_Alexa488),
    response.median = median(Intensity_MeanIntensity_Alexa488)
  ) ->
  KZ190n_pS1_T2_0.sum



cbind(KZ190n_pS1_T2_0.sum %>% dplyr::ungroup(),
      KZ190n_pS1_T2_0.sum %>% dplyr::filter(
        inhibitor.3.1 == 0,
        inhibitor.3.2 == 0) %>% 
        dplyr::ungroup() %>%
        dplyr::rename(scale.mean = response.mean, 
                      scale.median = response.median) %>%
        dplyr::select(scale.mean, scale.median) %>% 
        dplyr::ungroup() ) %>%
  dplyr::mutate(
    relative.mean = response.mean/scale.mean,
    relative.median = response.median/scale.median
  ) ->
  KZ190n_pS1_T2_0.sum

# color.limits <-
# quantile(x = KZ190n_pS1[["Intensity_MeanIntensity_Alexa488"]],
# probs = c(0.1,0.9))
# zakres: color.limits <- c(15, 300)
# mean-max: color.limits <- NULL

# color.limits <- c(15, max(KZ190n_pS1.sum$response.mean))
color.limits <- c(0,1)
g.inhibitors[["TYK2-ps1"]] <-
  ggplot(data = KZ190n_pS1_T2_0.sum,
         mapping = aes(x = factor(inhibitor.3.1),
                       y = factor(inhibitor.3.2),
                       fill = relative.median,
                       label = round(relative.median,2)
         )) +
  geom_tile(color = "white") +
  viridis::scale_fill_viridis(
    name = "Relative response\npY-STAT1",
    limits = color.limits) + 
  ggtitle("TYK2 = 0")+
  xlab("JAK1")+
  ylab("JAK2")+
  ITRC::theme_itrc() + 
  geom_text(size = 8)
g.inhibitors[["TYK2-ps1"]]



##### generic #####
pstat.type <- "pY-STAT1"
response.type <- "Intensity_MeanIntensity_Alexa488"
data <- KZ190n_pS1

pstat.type <- "pY-STAT3"
response.type <- "Intensity_MeanIntensity_Alexa555"
data <- KZ190n_pS3

scaled.by.inhib.val <- FALSE

data %>% 
  dplyr::group_by(
    inhibitor.3.1,
    inhibitor.3.2,
    inhibitor.3.3
  ) %>% 
  dplyr::summarise_(
    response.mean = paste("mean(", response.type,")" ),
    response.median = paste("median(", response.type,")" )
  ) ->
  data.sum

inhibitor.list <- c(
  "inhibitor.3.1",
  "inhibitor.3.2",
  "inhibitor.3.3")

inhibitors.labels <- c("JAK1", "JAK2", "TYK2")

g.inhibitor.list <- 
foreach(inhibitor.i = 1:length(inhibitor.list)) %do% {
  inhibitor_ = inhibitor.list[inhibitor.i]
  inhibitor.label <- inhibitors.labels[inhibitor.i]
  inhibitr.values.list <- (KZ190n_pS1 %>% dplyr::distinct_(inhibitor_) %>% dplyr::arrange_(inhibitor_))[[inhibitor_]]
  inhibitor.x <- inhibitor.list[inhibitor.list != inhibitor_][1]
  inhibitor.x.label <- inhibitors.labels[-inhibitor.i][1]
  inhibitor.y <- inhibitor.list[inhibitor.list != inhibitor_][2]
  inhibitor.y.label <- inhibitors.labels[-inhibitor.i][2]
  g.inhibitor.values.list <- 
    foreach(inhibitor.val_ = inhibitr.values.list) %do% {
      if(scaled.by.inhib.val){
        inhibitor.scale_ <- inhibitor.val_
      } else {
        inhibitor.scale_ <- 0
      }
      
      cbind(data.sum %>%
              dplyr::filter_(paste(inhibitor_, "==", inhibitor.val_)) %>%
              dplyr::ungroup(),
            data.sum %>%
              dplyr::filter_(paste(inhibitor_, "==", inhibitor.scale_)) %>% 
              dplyr::filter_(
                paste(inhibitor.x, "==", 0, "&",
                      inhibitor.y, "==", 0
                      ))%>% 
              dplyr::ungroup() %>%
              dplyr::rename(scale.mean = response.mean, 
                            scale.median = response.median) %>%
              dplyr::select(scale.mean, scale.median) %>% 
              dplyr::ungroup() ) %>%
        dplyr::mutate(
          relative.mean = response.mean/scale.mean,
          relative.median = response.median/scale.median
        ) ->
        data.sum.relative
      
      
  ggplot(data = data.sum.relative ,
       mapping = aes_string(x = paste("factor(", inhibitor.x, ")"),
                     y =  paste("factor(", inhibitor.y, ")"),
                     fill = "relative.median",
                     label = paste("round(relative.median,2)")
       )) +
  geom_tile(color = "white") +
  viridis::scale_fill_viridis(
    name = paste0("Relative response\n", pstat.type),
    limits = color.limits) + 
  ggtitle(paste(pstat.type, ":", inhibitor.label, "=", inhibitor.val_))+
  xlab(inhibitor.x.label)+
  ylab(inhibitor.y.label)+
  ITRC::theme_itrc() + 
  geom_text(size = 8)


    }
  g.inhibitor <- plot_grid(plotlist = g.inhibitor.values.list, nrow = 2, ncol =2 )
  return(g.inhibitor) 
}

g.inhibitor <- plot_grid(plotlist = g.inhibitor.list, ncol = 3)
ggsave(filename = paste(poster.path.list$output.fnp.dir, 
                        paste0("inhibitors-",
                               pstat.type
                               ,".pdf"), sep = "/"), 
       plot = g.inhibitor,
       width = 14*3,
       height = 12, useDingbats = FALSE
       )

#### Inhibition SCRC ####
bootstrap = TRUE
bootstrap.number = 32
bootstrap.sample_size = 1000
bootstrap.test.sample = TRUE
bootstrap.test.number <- 4
parallel_cores <- 8
KZ190n.zero <- subset(KZ190n, KZ190n$well.name %in% c( "D17", "D18"))

KZ190n_pS1 <- subset(KZ190n_filtered, KZ190n_filtered$antibody.1.1=="1:100")
pstat.type <- "pY-STAT1"
response.type <- "Intensity_MeanIntensity_Alexa488"
data <- KZ190n_pS1


pstat.type <- "pY-STAT3"
response.type <- "Intensity_MeanIntensity_Alexa555"
data <- KZ190n_pS3  


##

data %>%
  rbind(.,
        KZ190n_pS1.zero) ->
  data

inhibitor_ <- "inhibitor.3.3"

data %>% dplyr::filter(inhibitor.3.1 == 0, inhibitor.3.2 == 0) ->
  data.curve


data.curve %>%
  dplyr::distinct(stimulation.1.1, stimulation.1.2, inhibitor.3.3) ->
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
           dplyr::filter_("stimulation.1.1 != 0 & stimulation.1.2 != 0"))) ->
  data.curve.distinct
data.curve.distinct$level <- 1:nrow(data.curve.distinct)

data.curve %>%
  dplyr::left_join(data.curve.distinct,
                   by = c("stimulation.1.1", "stimulation.1.2", inhibitor_)) ->
  data.curve

levels.list <- 
  data.curve %>% 
  dplyr::filter(level != min(level),
                level != max(level)
                ) %>% 
  dplyr::arrange(level) %>%
  dplyr::distinct(level)
level.max <- max(data.curve$level)
level.min <- min(data.curve$level)

foreach(level_ = levels.list$level) %do% {
  
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
  data.frame(
    level = level_,
    complete = itrc.model$confusion.matrix[2,1], 
    partial = itrc.model$confusion.matrix[2,2]
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
    ) %>% dplyr::ungroup()
    
} %>% 
  do.call(what = rbind,
          args = .) %>%
  dplyr::left_join(data.curve.distinct,
                   by = "level") ->
  inhibitor.df


ggplot(data = inhibitor.df %>% dplyr::filter(type != "no_inhibition"),
       aes_string(
         x = inhibitor_,
         y = "prob",
         ymin = "prob - prob.sd",
         ymax = "prob + prob.sd",
         group = "type")) +
  geom_line() +
  geom_errorbar() +
  coord_cartesian(ylim = c(0,1)) +
  SysBioSigTheme::theme_sysbiosig()
  