### ###
### libraries graphics
### ###

filename.graphics       <- "libraries.R"
path.graphics  <- "R/graphics/"
files.graphics <- list.files(path = path.graphics,
                    pattern = ".R$", 
                    full.names = FALSE, 
                    recursive = TRUE)
files.graphics <- files.graphics[filename.graphics != files.graphics]

comment.graphics <- sapply(
  paste(path.graphics, files.graphics, sep = "/"),
  source,
  echo = FALSE,
  verbose = FALSE)

rm(filename.graphics,path.graphics,files.graphics,comment.graphics)