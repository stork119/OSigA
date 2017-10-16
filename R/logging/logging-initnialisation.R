### ###
### logging initialisation
### ###

InitLogging <- function(filename){
  print(filename)
  flog.appender(appender.file(filename))
  flog.threshold(INFO, name='logger.optimisation')
}