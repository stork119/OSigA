### ###
### logging initialisation
### ###

InitLogging <- function(filename){
  print(filename)
  flog.appender(appender.file(filename), name='logger.optimisation')
  flog.threshold(INFO, name='logger.optimisation')
}