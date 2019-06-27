### ###
### 2019-06-03-FNP
### ###

source("R/scripts/2017-10-05-KarolinaPoster/scriptsLibrary.R")
library(viridis)
library(SysBioSigTheme)
#### prepare data ####
poster.path.list <- getPathsList(normalisation = "ffc_filtered_oranized", date = '2018-03-29-experiments-poster_cyotkine')
poster.data.list <- readRDS(file = poster.path.list$rds.path)
poster.path.list$output.fnp.dir <- "resources/output/poster/2018-03-29-experiments-poster_cyotkine/FNPreport"
dir.create(path = poster.path.list$output.fnp.dir, recursive = TRUE)
#### pS1 priming vs non-primng ####
data.nonpriming <- poster.data.list$`2016-01-26-KA10-C` %>% dplyr::mutate(priming.1.2 = "U")
data.priming <- poster.data.list$`2016-01-26-KA10-IFNB`
ylim_ <- c(0,850)
data <- rbind(
  data.nonpriming,
  data.priming
)
stimulation.list <- "c(0, 0.05, 0.25, 5)"

#, "#a4d09e", "#7abd6f",  "#50af43", "#1fa637", "#177d34", "#175929"
colors.ifng.df <- data.frame(fill = c("#ffffff", "#bdd0ed", "#bdd0ed", "#7f9bd0", "#7e9fd3" ,"#6185c3", "#567abc", "#325aa6", "#213f72"),
           stimulation  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 5, 10),
           stringsAsFactors = FALSE) 
colors.ifng <- colors.ifng.df$fill  %>% as.vector()

colors.ifng_ifnb.df <- data.frame(fill = c("#c7e9c5", "#a4d09e", "#a4d09e", "#7abd6f",  "#50af43", "#1fa637", "#177d34", "#175929", "#0b2d14"),
                             stimulation  =  c(0, 0.01, 0.05, 0.1, 0.25,   0.5, 1, 5, 10),
                             stringsAsFactors = FALSE) 


names(colors.ifng) <- colors.ifng.df$stimulation
colors.ifnb.df <- data.frame(fill = c("#325aa6", "#175929"),
           priming.1.1  =  c(0, 1000),
           size  =  c(0.5, 1),
           stringsAsFactors = FALSE)  
colors.ifnb <- colors.ifnb.df$fill %>% as.vector()
names(colors.ifnb) <- colors.ifnb.df$priming.1.1 %>% as.vector()

colors.joined.df <-
  rbind(
    colors.ifng.df %>% dplyr::mutate(priming.1.1 = 0),
    colors.ifng_ifnb.df %>% dplyr::mutate(priming.1.1 = 1000)) %>% 
  dplyr::mutate(
    color_exp = paste(priming.1.1, stimulation)
  )
colors.joined <- colors.joined.df$fill
names(colors.joined) <- colors.joined.df$color_exp


column.names <- getColumnNames(data)

g.ps1 <- ggplot(data = 
         data %>% 
         dplyr::filter_(
           paste(column.names$stimulation, "%in%", stimulation.list, collapse = " ")) %>%
         dplyr::left_join(colors.ifng.df, by = "stimulation") %>%
         dplyr::left_join(colors.ifnb.df, by = "priming.1.1") %>%
          dplyr::left_join(colors.joined.df, by = c("stimulation", "priming.1.1")),
       mapping = aes_string(x = paste("factor(", column.names$stimulation, ")"),
                            y = column.names$response,
                            group = paste("interaction(", column.names$stimulation, ",", column.names$priming ,")"),
                            fill = paste("factor(", column.names$stimulation, ")"),
                            color = paste("factor(", column.names$priming, ")"))
) + 
  stat_boxplot(geom ='errorbar',   size = 0.75) + 
  geom_boxplot(outlier.size = 0.1,  size = 0.75, outlier.alpha = 0) +
  ITRC::theme_itrc() +
  coord_cartesian(ylim = ylim_) +
  scale_color_manual(
    name = "IFN-beta\nlevel [U]",
    values = colors.ifnb) +
  scale_fill_manual(
    name = "IFN-gamma\nlevel [ng/mL]",
    values = colors.ifng) +
  xlab("IFN-gamma level [ng/mL]") +
  ylab("cellular response")

ylim_ <- c(0,850)
g.ps1.diff_color <-
  ggplot(
    data = 
      data %>% 
      dplyr::filter_(
        paste(column.names$stimulation, "%in%", stimulation.list, collapse = " ")) %>%
      dplyr::left_join(colors.ifng.df, by = "stimulation") %>%
      dplyr::left_join(colors.ifnb.df, by = "priming.1.1") %>%
      dplyr::left_join(colors.joined.df, by = c("stimulation", "priming.1.1")),
    mapping = 
      aes_string(x = paste("factor(", column.names$stimulation, ")"),
                 y = column.names$response,
                 group = paste("interaction(", column.names$stimulation, ",", column.names$priming ,")"),
                 fill = "color_exp"
                 #, color = paste("factor(", column.names$priming, ")")
                 )) +
  stat_boxplot(geom ='errorbar',   size = 0.75) + 
  geom_boxplot(outlier.size = 0.1,  size = 0.75, outlier.alpha = 0) +
  ITRC::theme_itrc() +
  coord_cartesian(ylim = ylim_) +
  # scale_color_manual(
  #   name = "IFN-beta\nlevel [U]",
  #   values = colors.ifnb) +
  scale_fill_manual(
    name = "IFN-gamma\nlevel [ng/mL]",
    labels = c(0, 0.05, 0.25, 5, 0, 0.05, 0.25, 5),
    values = colors.joined) +
  xlab("IFN-gamma level [ng/mL]") +
  ylab("cellular response")

g.ps1.diff_color


ggsave(filename =paste(poster.path.list$output.fnp.dir, "IF-ps1.pdf", sep = "/"), 
       plot = g.ps1.diff_color,
       width = 8,
       height = 6, useDingbats = FALSE
)

#### SCRC pS1 priming vs non-primng ####
parallel_cores <- 8
bootstrap.sample_size = 1000
bootstrap.number = 8
bootstrap.test.sample <- NULL

model.nonpriming.itrc <- 
  ITRC::ITRC(data = data.nonpriming,
             signal = column.names$stimulation,
             response = column.names$response, 
             parallel_cores = parallel_cores,
             bootstrap.sample_size = bootstrap.sample_size,
             bootstrap.number = bootstrap.number,
             bootstrap.test.sample = bootstrap.test.sample
             )
model.priming.itrc <- 
  ITRC::ITRC(data = data.priming,
             signal = column.names$stimulation,
             response = column.names$response, 
             parallel_cores = parallel_cores,
             bootstrap.sample_size = bootstrap.sample_size,
             bootstrap.number = bootstrap.number,
             bootstrap.test.sample = bootstrap.test.sample
  )


g.itrc <- 
  ggplot(
    model.nonpriming.itrc$itrc %>% dplyr::mutate(priming.1.1 = 0) %>%
      rbind(.,
            model.priming.itrc$itrc %>% dplyr::mutate(priming.1.1 = 1000)),
    mapping = 
      aes_string(x = paste("factor(", column.names$stimulation, ")"),
                 y = "itrc", 
                 group = paste("factor(", column.names$priming, ")"),
                 color = paste("as.character(", column.names$priming, ")")
      )) + 
  geom_line(size = 1.5) +
  ITRC::theme_itrc() +
  coord_cartesian(ylim = c(0.5, 2.5)) + 
  scale_color_manual(values = colors.ifnb)


#### IRFS #### 
poster.path.list$fish.input <- "resources/input/poster/2018-03-29-experiments-poster_cyotkine/RNAFISH/"
experiments.list.fish <-
  data.frame(
    exp = c("2018-09-05-KZ-FISH12", "2018-10-18-KZ-FISH17", "2018-11-16-KZ-FISH27", "2018-11-29-KZ-FISH31", "2019-03-07-KZ-FISH39"),
    protein = c("Cxcl10", "Irf1", "Iigp1", "Parp14", "Tgtp2"), stringsAsFactors = FALSE
  )

foreach(exp.i = 1:nrow(experiments.list.fish)) %do% {
  path <- 
    paste(
      poster.path.list$fish.input,
      experiments.list.fish[exp.i,]$exp %>% as.character(),
      "CellsFiltered488.csv",
      sep = "/") 
  data.fish <- 
    (read.table(file = path, header = TRUE, sep = ",") %>% 
       SysBioSigTheme::normalize_data())$data
  data.fish %>% 
    dplyr::mutate(protein = as.character(experiments.list.fish[exp.i,"protein"])) %>%
    data.table::data.table() %>% 
      return()
} -> data.fish.list

names(data.fish.list) <- experiments.list.fish$protein

data.fish.list$Cxcl10 %>%
  dplyr::filter(time.1.1 == 120) ->
  data.fish.list$Cxxl10

data.fish.list$Parp14 %>%
  dplyr::filter(probe.1.1 == "Parp14") ->
  data.fish.list$Parp14

data.fish.list$Tgtp2 %>%
  dplyr::filter(probe.1.2 == "Tgtp2") ->
  data.fish.list$Tgtp2

data.fish.list$Iigp1 %>%
  dplyr::filter(probe.1.2 == "Iigp1") ->
  data.fish.list$Iigp1

foreach(data.fish = data.fish.list) %do% {
  data.fish %>%
  dplyr::select(
    protein,
    well.name,
    time.1.1,
    stimulation.1.1, stimulation.2.1, priming.1.1, priming.2.1,
    Intensity_MeanIntensity_Alexa546,
    Intensity_MedianIntensity_Alexa546
  ) %>% 
    data.table::data.table()
} -> 
  data.fish.list
data.fish.list %>%
  do.call(
    what = rbind,
    args = .) ->
   data.fish


data.fish %>% 
  dplyr::group_by(
    protein, time.1.1, stimulation.1.1, priming.1.1) %>%
  dplyr::summarise(
    intensity = median(Intensity_MeanIntensity_Alexa546)) %>%
  dplyr::filter(time.1.1 == 120) %>%
  dplyr::filter(stimulation.1.1 %in% c(0, 0.05, 0.25, 5)) ->
  data.fish.sum


# scale by priming 0 ifn 0 
data.fish.sum %>%
  dplyr::filter(stimulation.1.1 == 0, priming.1.1 == 0 ) %>%
  dplyr::group_by(
    protein) %>%
  dplyr::rename(scale = intensity) %>% 
  dplyr::ungroup() %>%
  dplyr::select(protein, scale) ->
  data.fish.sum.zero
  
data.fish.sum %>%
  dplyr::left_join(data.fish.sum.zero,
                   by = "protein") %>% 
  dplyr::mutate(response = intensity/scale) ->
  data.fish.sum
data.fish.sum %>% dplyr::ungroup() -> data.fish.sum

foreach(protein_ = 
          (data.fish.sum %>% 
             dplyr::distinct(protein))[["protein"]]) %do% {
  ggplot(
    data.fish.sum %>% dplyr::filter(protein == protein_),
    mapping = aes(
      x = factor(stimulation.1.1),
      y =  factor(priming.1.1),
      fill = response
    )
  ) +
    geom_tile() +
    ITRC::theme_itrc() +
    xlab("IFN-gamma level [ng/mL]") +
    ylab("IFN-beta\nlevel [U/mL]")
} ->
  g.list

names(g.list) <-
  (data.fish.sum %>% 
     dplyr::distinct(protein))[["protein"]]
# data.fish.sum %>% dplyr::filter(protein == "Irf1")

# g.list$Irf1 +     
#   scale_fill_viridis(
#     name = "cellular response rescaled", 
#     end =  0.8,
#     limits = c(0.9, 2.4)) 
#   
# g.list$Cxcl10 +     
#   scale_fill_viridis(
#     name = "cellular response rescaled", 
#     end =  0.8,
#     limits = c(0.9, 3.75)) 

library(cowplot)
plot_grid(plotlist = g.list, ncol = 1)
max(data.fish.sum$response)
g.list[[1]]


ggplot(
  data.fish.sum %>% dplyr::filter(protein != "Cxcl10"),
  mapping = aes(
    x = factor(stimulation.1.1),
    y =  interaction(-priming.1.1, protein),
    fill = response
  )
) +
  geom_tile() +
  ITRC::theme_itrc() +
  scale_fill_viridis()


foreach(protein_ = 
          (data.fish.sum %>% 
             dplyr::distinct(protein))[["protein"]]) %do% {
               ggplot(
                 data.fish.sum %>% dplyr::filter(protein == protein_),
                 mapping = aes(
                   x = factor(stimulation.1.1),
                   y =  response,
                   group = priming.1.1,
                   color = factor(priming.1.1)
                 )
               ) +
                 geom_line(size = 1) +
                 ITRC::theme_itrc() +
                 xlab("IFN-gamma level [ng/mL]") +
                 ylab("relative response")
             } ->
  g.line.list

names(g.line.list) <-
  (data.fish.sum %>% 
     dplyr::distinct(protein))[["protein"]]


ggsave(filename =paste(poster.path.list$output.fnp.dir, "RNAFISH-tgtp2.pdf", sep = "/"), 
       plot = g.line.list$Tgtp2  +
         coord_cartesian(ylim = c(1,1.5)) + 
         scale_color_manual(values = colors.ifnb,
                            name = "IFN-beta\nlevel [U/mL]"),
       width = 6,
       height = 3, useDingbats = FALSE
)

ggsave(filename =paste(poster.path.list$output.fnp.dir, "RNAFISH-parp14.pdf", sep = "/"), 
       plot = g.line.list$Parp14  +
         coord_cartesian(ylim = c(1,1.5)) + 
         scale_color_manual(values = colors.ifnb,
                            name = "IFN-beta\nlevel [U/mL]"),
       width = 6,
       height = 3, useDingbats = FALSE
)

ggsave(filename =paste(poster.path.list$output.fnp.dir, "RNAFISH-Iigp1.pdf", sep = "/"), 
       plot = g.line.list$Iigp1  +
         coord_cartesian(ylim = c(1,1.5)) + 
         scale_color_manual(values = colors.ifnb,
                            name = "IFN-beta\nlevel [U/mL]"),
       width = 6,
       height = 3, useDingbats = FALSE
)


ggsave(filename =paste(poster.path.list$output.fnp.dir, "RNAFISH-irf1.pdf", sep = "/"), 
       plot = g.line.list$Irf1  +
         coord_cartesian(ylim = c(1,2.5)) + 
         scale_color_manual(values = colors.ifnb,
                            name = "IFN-beta\nlevel [U/mL]"),
       width = 6,
       height = 3, useDingbats = FALSE
)


ggsave(filename =paste(poster.path.list$output.fnp.dir, "RNAFISH-cxcl10.pdf", sep = "/"), 
       plot = g.line.list$Cxcl10 +
         coord_cartesian(ylim = c(1,3.5)) + 
         scale_color_manual(values = colors.ifnb,
                            name = "IFN-beta\nlevel [U/mL]"),
       width = 6,
       height = 3, useDingbats = FALSE
)
#### SCRC downstreams ####
bootstrap = TRUE
bootstrap.number =  64
bootstrap.sample_size = 1000
bootstrap.test.sample = TRUE
bootstrap.test.number = 4
parallel_cores = 8

stimulation.list <- c(0,0.05, 0.25, 5)
data.fish %>% 
  dplyr::filter(stimulation.1.1 %in% stimulation.list) ->
  data.fish

title_ = ""
xlab_ = ""
ylab_ = ""
ylimits_ <- c(0,3)

data.fish %>% dplyr::distinct(protein)
protein_ <- "Tgtp2"
response_ <- "Intensity_MeanIntensity_Alexa546"
signal_ <- "stimulation.1.1"
scrc.model <- list()
scrc.plot.list  <- list()
scrc.plot  <- list()
for(protein_ in (data.fish %>% dplyr::distinct(protein))[["protein"]]) {
  scrc.model[[protein_]] <- 
  list( 
    nopriming = 
      ITRC::ITRC(
        data = data.fish %>% 
          dplyr::filter(protein == protein_) %>% 
          dplyr::filter(priming.1.1 == 0), 
        signal = signal_,
        response = response_,
        bootstrap.number = bootstrap.number,
        bootstrap = bootstrap,
        bootstrap.sample_size = bootstrap.sample_size,
        bootstrap.test.sample = bootstrap.test.sample,
        bootstrap.test.number = bootstrap.test.number,
        parallel_cores = parallel_cores
      ),
    priming = 
      ITRC::ITRC(
        data = data.fish %>% 
          dplyr::filter(protein == protein_) %>% 
          dplyr::filter(priming.1.1 == 1000), 
        signal = signal_,
        response = response_,
        bootstrap.number = bootstrap.number,
        bootstrap = bootstrap,
        bootstrap.sample_size = bootstrap.sample_size,
        bootstrap.test.sample = bootstrap.test.sample,
        bootstrap.test.number = bootstrap.test.number,
        parallel_cores = parallel_cores
      ))
}

theme.signal <-
  ITRC::GetRescaledSignalTheme(
    model = scrc.model$Irf1$nopriming,
    rescale.fun = function(x){log(x = x, base = 10)},
    pallete.args = 
      list(end = 0.95))

for(protein_ in (data.fish %>% dplyr::distinct(protein))[["protein"]]) {
scrc.plot.list[[protein_]] <-
  list(
    nopriming = 
      ITRC::plotITRCWaves(
        model = scrc.model[[protein_]]$nopriming, 
        ylimits_ = ylimits_,
        theme.signal = theme.signal,
        title_ = paste(protein_, "no-priming"),
        xlab_ = xlab_,
        ylab_ = ylab_),
    priming = 
      ITRC::plotITRCWaves(
        model = scrc.model[[protein_]]$priming, 
        ylimits_ = ylimits_,
        theme.signal = theme.signal,
        title_ = paste(protein_, "priming"),
        xlab_ = xlab_,
        ylab_ = ylab_)
  )
  scrc.plot[[protein_]] <-
    plot_grid(plotlist = scrc.plot.list[[protein_]], nrow = 1)
}

scrc.plot.combined <- 
  plot_grid(plotlist = scrc.plot,
            ncol = 1)


saveRDS(object = scrc.model,
        file = paste(poster.path.list$output.fnp.dir, "RNAFISH-SCRC.rds", sep = "/"))

ggsave(filename =paste(poster.path.list$output.fnp.dir, "RNAFISH-SCRC-plot.pdf", sep = "/"), 
       plot = scrc.plot.combined,
       width = 12,
       height = 30, useDingbats = FALSE
)
