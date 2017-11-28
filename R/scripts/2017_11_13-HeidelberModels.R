library(ggplot2)
n <- 1000
parameters <- rnorm(n=n)
parameters.df <- data.frame(par = parameters)

ggplot(parameters.df, mapping = aes(x = par)) + geom_density()


par <- 0.5
data.df <- data.frame(Y = rnorm(n = n, mean =  par), par = par)
par <- -2
data.df <- rbind(
  data.df,
  data.frame(Y = rnorm(n = n, mean =  par),
             par = par)
)

ggplot(data.df, mapping = aes(x = Y, group = par)) + geom_density()


ggplot(data.df, mapping = aes(x = Y)) + geom_density()