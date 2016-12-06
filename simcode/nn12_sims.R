library(pevosim)
library(dplyr)
library(tidyr)

pardf <- expand.grid(dt=c(1,1/10,1/20,1/50),
                     R0=c(10,1e6),
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
    save("pardf","resList",file="nn12.rda")
}
{
if (FALSE)  {
library(ggplot2); theme_set(theme_bw())
rL <- lme4:::namedList(base,lo_dt,base_hi_dt,tstep_hi_dt)
x1 <- rL %>% bind_rows(.id="type") %>%
    select(type,time,mean_R0,sd_R0) %>%
    gather(var,val,-c(type,time))

ggplot(x1,aes(time,val,colour=type))+
    geom_line()+facet_wrap(~var,scale="free")+
    scale_y_log10()
}
