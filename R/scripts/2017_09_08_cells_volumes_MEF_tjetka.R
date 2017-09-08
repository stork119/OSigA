### ###
### 2017_09_08_cells_volumes_MEF_tjetka
### sources : resources/input/cells_volumes/2017-09-08-MEF_tjetka.txt
### ###

source("R/libraries.R")

# data.list <- list()
data.list$cellsvolumes <- read.table(file = "resources/input/cells_volumes/2017-09-08-MEF_tjetka.txt",
                                     header = TRUE)

GetAreaFromVolume <- function(
  V,
  S0 = (36*pi)^(1/3)){
  return(S0*(V^(2/3)))
}


data.list$cellsvolumes <- 
  data.list$cellsvolumes %>%
  dplyr::mutate(CellsArea = GetAreaFromVolume(Cells),
                NucleiArea = GetAreaFromVolume(Nuclei))

write.table(x = data.list$cellsvolumes,
            file = "resources/input/cells_volumes/2017-09-08-MEF_tjetka.csv",
            sep = ",",
            col.names = TRUE,
            row.names = FALSE)

data.list$cellsvolumes.summarise <-
  data.list$cellsvolumes %>%
  dplyr::summarise(Cells = mean(Cells),
                   CellsArea = mean(CellsArea),
                   Nuclei = mean(Nuclei),
                   NucleiArea = mean(NucleiArea))

  