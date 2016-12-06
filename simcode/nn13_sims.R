library(pevosim)
library(dplyr)
library(tidyr)

pardf13 <- expand.grid(dt=c(1/100,1/200,1/500),
                     R0=c(10,1e7),
                     N=c(30,100,300,1000))

resList13 <- vector(nrow(pardf13),mode="list")
for (i in 1:nrow(pardf13)) {
    print(pardf13[i,])
    resList13[[i]] <-
        try(with(as.list(pardf13[i,]),
             run_sim(nt=1e5,seed=101,dt=dt,
                     mod_init=list(R0=R0),
                     mod_params=list(N=N),
                     progress=TRUE,
                     hazard=TRUE)))
    save("pardf13","resList13",file="nn13.rda")
}
