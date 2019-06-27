data.nonpriming <- poster.data.list$`2016-01-26-KA10-C` %>% dplyr::mutate(priming.1.2 = "U") 
data.priming <- poster.data.list$`2016-01-26-KA10-IFNB`

data.nonpriming <- poster.data.list$`2016-01-28-KA11-C` %>% dplyr::mutate(priming.1.2 = "U")
data.priming <- poster.data.list$`2016-01-28-KA11-IFNB`

data <- rbind(
  data.nonpriming,
  data.priming
) %>% 
  dplyr::filter(time.1.1 == 15) %>%
  data.table()  %>%
  rbind(.,
        (data.nonpriming %>%
           dplyr::filter(time.1.1 == 0, stimulation.1.1 == 0.1) %>% 
           dplyr::mutate(stimulation.1.1 = 0, time.1.1 = 15, stimulation = 0)))
data %>% 
  dplyr::mutate_(response.log = paste("log(", column.names$response, ")")) ->
  data

#xlim_ <- c(0,600)
fill.values <- c("#325aa6",  "#175929")
names(fill.values) <- c(0,1000)

g.histogram <-
  ggplot(data = 
    data %>% dplyr::filter(stimulation %in% c(0,0.05, 0.25, 5)),
         # dplyr::filter_(
         #   paste(column.names$stimulation, "%in%", stimulation.list, collapse = " ")) %>%
         # dplyr::left_join(colors.ifng.df, by = "stimulation") %>%
         # dplyr::left_join(colors.ifnb.df, by = "priming.1.1") %>%
         # dplyr::left_join(colors.joined.df, by = c("stimulation", "priming.1.1")),
       mapping = aes_string(#x = paste("factor(", column.names$stimulation, ")"),
                            x =  "log(Intensity_MeanIntensity_Alexa)", 
                            #group = paste("interaction(", column.names$stimulation, ",", column.names$priming ,")"),
                            # fill = paste("factor(", column.names$stimulation, ")"),
                            fill = paste("factor(", column.names$priming, ")")
       )
) + 
  geom_histogram(
    binwidth = 0.025,
    position = 
      "dodge",
    alpha = 0.75
    ) +
  SysBioSigTheme::theme_sysbiosig() +
  facet_grid(stimulation.1.1~.) +
  scale_fill_manual(values = fill.values)


ggsave(filename = paste(poster.path.list$output.fnp.dir, 
                        paste0("ifn-histograms.pdf"), sep = "/"), 
       plot = g.histogram,
       width = 12,
       height = 8, useDingbats = FALSE
)

#####$
#devtools::install_github("stork119/SysBioSigTheme", auth_token =  "46a792807a855397b2fc325919db10179050c7f8")
library(SysBioSigTheme)
library(data.table)
# data.nonpriming <- poster.data.list$`2016-01-26-KA10-C` %>% dplyr::mutate(priming.1.2 = "U") 
# data.priming <- poster.data.list$`2016-01-26-KA10-IFNB`

data.nonpriming <- poster.data.list$`2016-01-28-KA11-C` %>% dplyr::mutate(priming.1.2 = "U")
data.priming <- poster.data.list$`2016-01-28-KA11-IFNB`

data <- rbind(
  data.nonpriming,
  data.priming
) %>% 
  data.table() 


g.scatterboxplot <-
  SysBioSigTheme::ScatterBoxplotGGplot(
    data = data  %>% 
      dplyr::mutate_(response.log = paste("log(", column.names$response, ")")) %>% 
      dplyr::filter(stimulation.1.1 == 5),
    x_ = "time.1.1",
    y_ = "response.log",
    facet.cols = "priming.1.1"
    )

ggsave(filename = paste(poster.path.list$output.fnp.dir, 
                        paste0("ifn-scaterboxplot-time.pdf"), sep = "/"), 
       plot = g.scatterboxplot,
       width = 12,
       height = 6, useDingbats = FALSE
)

#### SCRC ####

data <- rbind(
  data.nonpriming,
  data.priming
) %>% 
  dplyr::filter(time.1.1 == 15) %>%
  data.table()  %>%
  rbind(.,
        (data.nonpriming %>%
           dplyr::filter(time.1.1 == 0, stimulation.1.1 == 0.1) %>% 
           dplyr::mutate(stimulation.1.1 = 0, time.1.1 = 15, stimulation = 0)))
data %>% 
  dplyr::mutate_(response.log = paste("log(", column.names$response, ")")) ->
  data

bootstrap = TRUE
bootstrap.number = 64
bootstrap.sample_size = 1000
bootstrap.test.sample = TRUE
bootstrap.test.number <- 4
parallel_cores <- 8
xlab_ <- ""
ylab_ <- ""
title_ <- ""
ylimits_ <- c(0,4)
itrc.nonpriming <- ITRC::ITRC(
  data = data %>% dplyr::filter(priming.1.1 == 0), 
  signal = column.names$stimulation,
  response = column.names$response,
  bootstrap.number = bootstrap.number,
  bootstrap = bootstrap,
  bootstrap.sample_size = bootstrap.sample_size,
  bootstrap.test.sample = bootstrap.test.sample,
  bootstrap.test.number = bootstrap.test.number,
  parallel_cores = parallel_cores
)

itrc.priming <- ITRC::ITRC(
  data = data %>% dplyr::filter(priming.1.1 == 1000), 
  signal = column.names$stimulation,
  response = column.names$response,
  bootstrap.number = bootstrap.number,
  bootstrap = bootstrap,
  bootstrap.sample_size = bootstrap.sample_size,
  bootstrap.test.sample = bootstrap.test.sample,
  bootstrap.test.number = bootstrap.test.number,
  parallel_cores = parallel_cores
)

theme.signal <-
  ITRC::GetRescaledSignalTheme(
    model = itrc.nonpriming,
    rescale.fun = function(x){log(x = x, base = 10)},
    pallete.args = 
      list(end = 0.95))
g.ifn.scrc.list<-list()
g.ifn.scrc.list[["nonpriming"]] <- 
  ITRC::plotITRCWaves(
  model = itrc.nonpriming, 
  ylimits_ = ylimits_,
  theme.signal = theme.signal,
  title_ = title_,
  xlab_ = xlab_,
  ylab_ = ylab_)

g.ifn.scrc.list[["priming"]] <-
ITRC::plotITRCWaves(
  model = itrc.priming,
 ylimits_ = ylimits_,
  theme.signal = theme.signal,
title_ = title_,
xlab_ = xlab_,
ylab_ = ylab_)

g.ifn.scrc  <- cowplot::plot_grid(plotlist = g.ifn.scrc.list, nrow = 1)
ggsave(filename = paste(poster.path.list$output.fnp.dir, 
                        paste0("ifn-SCRC.pdf"), sep = "/"), 
       plot = g.ifn.scrc,
       width = 12,
       height = 6, useDingbats = FALSE
)