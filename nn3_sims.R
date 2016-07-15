## nearly neutral sims ... Latin hypercube

source("nn2_funs.R")
do_runs <- FALSE
nsim <- 1000
lhs_df <- cbind(
    mu=10^seq(-3,-1,length.out=nsim),
    mut_mean=seq(-4,-0.5,length.out=nsim),
    mut_sd=10^seq(-1,0.5,length.out=nsim),
    N=round(10^seq(0.5,3,length.out=nsim)),
    gamma=1/round(seq(2,10,length.out=nsim)))
set.seed(101)
for (i in 2:ncol(lhs_df)) {
    lhs_df[,i] <- sample(lhs_df[,i])
}
## pairs(lhs_df,pch=".",gap=0)
fn <- "nn_runs_3.rda"
if (do_runs) {
 resList <- list()
 for (i in 1:nrow(lhs_df)) {
     pars <- lhs_df[i,]
     cat(c(i,pars,"\n")
         argList <- c(list(nt=1e6,rptfreq=1000,R0_init=4,seed=101,
                           progress=TRUE),
                      as.list(pars))
         resList[[i]] <- try(do.call(run_sim,argList))
         save("resList",file=fn)
     }
 }
} else {
    load(fn)
    length(resList)
    library(plyr)
    library(dplyr)
    library(tidyr)
    rc1 <- ldply(resList,
          function(x) {
        x %>% filter(time>9e5) %>%
            summarise(
                I=mean(I),
                mean_lbeta=mean(mean_lbeta),
                sd_lbeta=mean(sd_lbeta))
        })
    rc2 <- cbind(lhs_df[1:nrow(rc1),],rc1)
    library(ggplot2); theme_set(theme_bw()); library(viridis)
    rc3 <- rc2 %>% mutate(Iprop=I/N,run=1:nrow(rc2)) %>% select(-I)
    rc4A <- rc3 %>% select(run,mu,mut_mean,mut_sd,N,gamma) %>%
        gather(var1,val1,-run)
    rc4B <- rc3 %>% select(-mu,-mut_mean,-mut_sd,-N,-gamma) %>%
        gather(var2,val2,-run)
    lhs_df <- data.frame(run=1:nrow(lhs_df),lhs_df)
    rc5 <- full_join(rc4A,rc4B,by="run")
    rc6 <- rc5 %>% full_join(lhs_df) %>% na.omit()
    ggplot(rc6,aes(val1,val2))+
        geom_point(aes(colour=mut_mean))+
        facet_grid(var2~var1,scale="free")+    
        scale_color_viridis()
    library(GGally)
    ggpfun <- function(var="mean_lbeta") {
        ggp1 <- ggpairs(rc3,
                        mapping = ggplot2::aes_string(colour = var),
                        columns=1:5,
                        lower = list(continuous = wrap("points",size=1)), ## alpha = 0.3,size=0.5)),
                        diag = list(continuous = "blankDiag"),
                        upper = list(continuous = "blank"))
        return(trim_gg(tweak_colours_gg(ggp1)))
    }
    source("ggpairsfuns.R")

    gList <- lapply(c("mean_lbeta","Iprop","sd_lbeta"),ggpfun)
    print(gList[[3]],spacingProportion=0)

