library(adagio)

par <- par.def
par.lower <- par*10^(-2)
par.upper <- par*10^(2)
maxit <- 1000

path.cmaes <- paste(path.output, "cmaes", sep = "")
dir.create(path.cmaes, recursive = TRUE)

cma_es.res <- pureCMAES(par = par, 
                   fun = optimisation,
                   lower = par.lower,
                   upper = par.upper,
                    stopeval = maxit,
                   #control = list(maxit = maxit, max.call = maxit), 
                   fun_run_model = run_model_mean,
                   variables = variables,
                   tmesh = tmesh,
                   tmesh.list = tmesh.list,
                   stimulation.list = stimulation.list,
                   background = background,
                   data.exp.grouped = data.exp.grouped,
                   fun.likelihood = fun.likelihood.mvn.mean)  

par <- scan('')
path.cmaes.opt <- paste(path.cmaes, "2017-01-26", sep ="/")

#### ####
res.par$data.model$likelihood <- likelihood(data.model = res.par$data.model,
                                            data.exp.grouped = data.exp.grouped,
                                            fun.likelihood = fun.likelihood.mvn)

res.def$data.model$likelihood <- likelihood(data.model = res.def$data.model,
                                            data.exp.grouped = data.exp.grouped,
                                            fun.likelihood = fun.likelihood.mvn)

data.model <- res.def$data.model
data.model <- res.par$data.model
data.model.i <- 49

sapply(1:nrow(data.model),
       function(data.model.i){
         data.model.tmp <- data.model[data.model.i,]
         return((
            d.removed %>%
         mutate(likelihood = 
                            do.call(fun.likelihood,list(logintensity, 
                                                        intensity, 
                                                        data.model.tmp = data.model.tmp))) %>%
                   summarise(likelihood.sum = 
                               sum(likelihood))$likelihood.sum)
       }
)

         
         d <- data.exp.grouped %>%
           filter(priming == data.model.tmp$priming,
                  time == data.model.tmp$time,
                  stimulation == data.model.tmp$stimulation)
         
         d.quantile <- quantile(d$intensity)
         d.quantile[2]
         
         d.removed <- d %>%
           filter(intensity > d.quantile[2],
                  intensity < d.quantile[4])
         