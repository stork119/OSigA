data.model <- dm$data.model

plot.args$theme.title_size <-36
plot.args$legend.position <-"center"
gplot <- ggplot(data.model %>%  
                  dplyr::mutate(type = "model") %>% 
                  dplyr::filter(stimulation %in% c(0.01, 0.25, 5)),
                mapping = aes(x = time, y = log(m.norm), group = type, color = type)) +
  # geom_point() +
  
  facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
  # geom_boxplot(data = data.list$data.exp.norm  %>% 
  #                dplyr::filter(stimulation %in% c(0.01, 0.25, 5)) %>% 
  #                dplyr::mutate(type = "data"),
  #              mapping = aes(x = time, y = logintensity, group = interaction(priming, time, stimulation)),
  #              color = "black") +
  geom_line(data = data.exp.summarise.optimisation %>% 
              mutate(type = "data")%>% 
              dplyr::filter(stimulation %in% c(0.01, 0.25, 5)),
            color = "black") +
  geom_errorbar(data = data.exp.summarise.optimisation %>%
                  mutate(type = "data")%>% 
                  dplyr::filter(stimulation %in% c(0.01, 0.25, 5)),
                mapping = aes(ymin = mean.lmvn - sqrt(sd.lmvn),
                              ymax = mean.lmvn + sqrt(sd.lmvn)),
                color = "black") +
  geom_line(size = 1.25) +
  ylab("log(cellular response)") + 
  xlab("time [min]") + 
  ggtitle("Cells MEF. Cellular response in IFN-gamma")


plot.args.ggsave$width <- 16
do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.list$optimisation.analysis, "best-model-3", "models_compare_GM.pdf", sep = "/"),
                           plot = gplot)))

gplot <- ggplot(data.model %>%  
                  dplyr::mutate(type = "model") %>% 
                  dplyr::filter(stimulation %in% c(0.01, 0.25, 5)),
                mapping = aes(x = time, y = log(m.norm), group = type, color = type)) +
  # geom_point() +
  
  facet_grid(priming ~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
  # geom_boxplot(data = data.list$data.exp.norm  %>% 
  #                dplyr::filter(stimulation %in% c(0.01, 0.25, 5)) %>% 
  #                dplyr::mutate(type = "data"),
  #              mapping = aes(x = time, y = logintensity, group = interaction(priming, time, stimulation)),
  #              color = "black") +
  geom_line(data = data.exp.summarise.optimisation %>% 
              mutate(type = "data")%>% 
              dplyr::filter(stimulation %in% c(0.01, 0.25, 5)),
            color = "black") +
  geom_errorbar(data = data.exp.summarise.optimisation %>%
                  mutate(type = "data")%>% 
                  dplyr::filter(stimulation %in% c(0.01, 0.25, 5)),
                mapping = aes(ymin = mean.lmvn - sqrt(sd.lmvn),
                              ymax = mean.lmvn + sqrt(sd.lmvn)),
                color = "black") +
  geom_line(size = 1.25) +
  ylab("log(cellular response)") + 
  xlab("time [min]") + 
  ggtitle("Cells MEF. Cellular response in IFN-gamma")


plot.args.ggsave$width <- 16
do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.list$optimisation.analysis, "best-model-3", "models_compare_errorbar_GM.pdf", sep = "/"),
                           plot = gplot)))

gplot.data <- ggplot(data = data.list$data.exp.norm  %>% 
                       dplyr::filter(stimulation %in% c(1)) %>%
                       dplyr::filter(priming %in% c(0)) %>% 
                       dplyr::mutate(type = "data"),
                     mapping = aes(x = time, y = logintensity, group = position),#interaction(priming, time, stimulation)),
                     color = "black") +
  geom_boxplot() +
  facet_grid(~ stimulation) + do.call(what = theme_jetka, args = plot.args) +
  ylab("log(cellular response)") + 
  xlab("time [min]") + 
  ggtitle("Cells MEF. Cellular response on 1 ng/ml of IFN-gamma")


plot.args.ggsave$width <- 16
plot.args.ggsave$height <-8
do.call(what = ggsave,
        args = append(plot.args.ggsave,
                      list(filename = paste(path.list$optimisation.analysis, "best-model-3", "data_GM.pdf", sep = "/"),
                           plot = gplot.data)))

