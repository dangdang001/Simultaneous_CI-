install.packages('doParallel')
library(doParallel)
system.time(foreach(i=1:10000) %do% sum(tanh(1:i)))


registerDoParallel(cores=10)

system.time(foreach(i=1:10000) %dopar% sum(tanh(1:i)))
