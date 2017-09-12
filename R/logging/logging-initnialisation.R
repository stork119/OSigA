### ###
### logging initialisation
### ###

#logReset()
# basicConfig(level='DEBUG')
# addHandler(writeToFile, file="~/Documents/modelling/scripts/karol.log", logger='logger.optimisation', level = 'DEBUG')
# #setLevel('DEBUG',  getHandler('writeToFile'))
# setLevel('DEBUG',  getHandler('logger.optimisation'))

InitLogging <- function(filename){
  paste(filename)
  flog.appender(appender.file(filename))
  flog.threshold(INFO, name='logger.optimisation')
}