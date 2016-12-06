library(pevosim)
library(dplyr)
library(tidyr)

pardf <- expand.grid(dt=c(1/100,1/200,1/500),
                     R0=c(10,1e7),
                     N=c(30,100,300,1000))

resList <- vector(nrow(pardf),mode="list")
for (i in 1:nrow(pardf)) {
    print(pardf[i,])
    resList[[i]] <-
        try(with(as.list(pardf[i,]),
             run_sim(nt=1e5,seed=101,dt=dt,
                     mod_init=list(R0=R0),
                     mod_params=list(N=N),
                     progress=TRUE,
                     hazard=TRUE)))
    save("pardf13","resList13",file="nn13.rda")
}
