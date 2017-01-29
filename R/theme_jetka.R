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
theme_jetka <-  function (base_size = 12,
                          base_family = "sans",
                          text_size = 18,
                          text_size_2 = 3*text_size/4) 
{
  (theme_foundation(base_size = base_size, base_family = base_family) + 
     theme(line = element_line(), 
           rect = element_rect(fill = ggthemes_data$fivethirtyeight["ltgray"], linetype = 0, colour = NA),
           plot.title = element_text(colour = ggthemes_data$fivethirtyeight["dkgray"], vjust = 1, hjust = 0.5, size=text_size, face="bold"), 
           text = element_text(colour = ggthemes_data$fivethirtyeight["dkgray"]), 
           axis.title = element_text(size=text_size),
           axis.title.y = element_text(angle=90),
           axis.text = element_text(size=text_size_2, face="bold"),
           axis.text.x = element_text(angle = 90),
           axis.ticks = element_blank(),
           axis.line.x = element_line(),
           axis.line.y = element_blank(), 
           legend.position="right",
           legend.background = element_rect(fill = "white"),
           panel.grid = element_line(colour = NULL), 
           panel.grid.major = element_line(colour = ggthemes_data$fivethirtyeight["medgray"]), 
           panel.grid.minor = element_blank(), 
           panel.background = element_rect(fill = "white"),
           # plot.title = element_text(hjust = 0, size = rel(1.75), face = "bold"), 
           plot.margin = unit(c(1,1, 1, 1), "lines"), 
           plot.background = element_rect(fill = "white"),
           
           strip.background = element_rect(fill = "white"),
           strip.text = element_text(size= text_size_2, face="bold", vjust = 0.5, lineheight = text_size_2*3)))
}
