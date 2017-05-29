### ###
###
###
###
### ###

#### PREPROCESSING  ####

# ### libraries ###
# try({package.list <- list("ggplot2", "ggthemes", "grid")
#   package.load <- sapply(package.list, function(package.name){
#   package.exist <- require(package.name, character.only = TRUE)
#   if(!package.exist){
#     install.packages(package.name)
#     return(library(package.name, character.only = TRUE))
#   }
#   return(package.exist)
# })
# })
# 
# ### sources ###
# # wd.tmp <- "" ### Rstudio 
# wd.tmp <- dirname(sys.frame(1)$ofile) ### script

#### MAIN  ####
theme_jetka <-  function (theme.base_size = 12,
                          theme.base_family = "sans",
                          theme.title_size = 18,
                          theme.text_size = 3*theme.title_size/4,
                          theme.margins = c(1,1, 1, 1),
                          ...) 
{
  (theme_foundation(base_size = theme.base_size, base_family = theme.base_family) + 
     theme(line = element_line(), 
           rect = element_rect(fill = 
                                 ggthemes_data$fivethirtyeight["ltgray"], 
                               linetype = 0, colour = NA),
           plot.title = element_text(colour = 
                                       ggthemes_data$fivethirtyeight["dkgray"],
                                     vjust = 1, hjust = 0.5, 
                                     size=theme.title_size,
                                     face="bold"), 
           text = element_text(colour = ggthemes_data$fivethirtyeight["dkgray"]), 
           axis.title = element_text(size=theme.title_size),
           axis.title.y = element_text(angle=90),
           axis.text = element_text(size=theme.text_size, face="bold"),
           axis.text.x = element_text(angle = 90),
           axis.ticks = element_blank(),
           axis.line.x = element_line(),
           axis.line.y = element_blank(), 
           legend.position="right",
           legend.background = element_rect(fill = "white"),
           panel.grid = element_line(colour = NULL), 
           panel.grid.major = 
             element_line(colour = ggthemes_data$fivethirtyeight["medgray"]), 
           panel.grid.minor = 
             element_blank(), 
           panel.background = element_rect(fill = "white"),
           # plot.title = element_text(hjust = 0, size = rel(1.75), face = "bold"), 
           plot.margin = unit(theme.margins, "lines"), 
           plot.background = element_rect(fill = "white"),
           
           strip.background = element_rect(fill = "white"),
           strip.text = 
             element_text(size= theme.text_size,
                          face="bold", 
                          vjust = 0.5, 
                          lineheight = theme.text_size*3)))
}
