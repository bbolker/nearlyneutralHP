library(dplyr)
load("../simdata/nn12.rda")
load("../simdata/nn13.rda")
cList <- c(resList,resList13)
cdf <- rbind(pardf,pardf13)
nrun <- length(cList)
pardf2 <- cdf %>%
    mutate(run=as.character(1:nrun))
names(cList) <- 1:nrun
library(ggplot2); theme_set(theme_bw())

x1 <- cList %>% bind_rows(.id="run") %>%
    select(run,time,mean_R0) %>%
    full_join(pardf2,by="run") %>%
    rename(start_R0=R0) %>%
    na.omit %>%
    gather(var,val,-c(run,dt,start_R0,N,time))

## all runs equilibrate (more or less)
ggplot(x1,aes(time/1000,val,colour=factor(start_R0)))+
    geom_line()+facet_grid(N~dt,scale="free",
                           labeller=label_both)+
    scale_colour_brewer(palette="Set1")+
    theme(panel.spacing=grid::unit(0,"lines"))+
    scale_y_log10()

x1s <- x1 %>% filter(time>8e4) %>%
    group_by(dt,start_R0,N) %>%
    summarise(mean_R0=mean(val),
              sd_R0=sd(val))

## for a given dt, logarithmic R0 increase with N
ggplot(x1s,aes(N,mean_R0,
               colour=factor(start_R0),
               shape=factor(dt)))+
    geom_pointrange(aes(ymin=mean_R0-2*sd_R0,
                        ymax=mean_R0+2*sd_R0))+
    geom_line(aes(group=interaction(factor(dt),factor(start_R0))))+
    scale_x_log10()+
    scale_y_log10()

## for a given dt, logarithmic R0 increase with 1/dt ... ??
ggplot(x1s,aes(1/dt,mean_R0,
               colour=factor(N),
               shape=factor(start_R0)))+
    geom_pointrange(aes(ymin=mean_R0-2*sd_R0,
                        ymax=mean_R0+2*sd_R0))+
    geom_line(aes(group=interaction(factor(N),factor(start_R0))))+
    scale_y_log10()+
    scale_x_log10()



