### ###
### rbind.smart
### ###


rbind.smart <- function(
  data.1, 
  data.2,
  ...){
  data.1.c <- colnames(data.1)
  data.2.c <- colnames(data.2)
  cols <- intersect(data.1.c, data.2.c)
  rbind(data.1[,cols], data.2[,cols], ...)
}