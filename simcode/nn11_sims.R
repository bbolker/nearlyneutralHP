library(pevosim)

t1A <- system.time(res1A <- run_sim(nt=1e4,rptfreq=10,seed=101,progress=TRUE))
t1B <- system.time(res1B <- run_sim(nt=1e4,rptfreq=10,seed=101,discrete=FALSE,
                                    progress=TRUE))
t1C <- system.time(res1C <- run_sim(nt=1e4,rptfreq=10,seed=101,progress=TRUE,
                           dt=1,hazard=TRUE))

t1D <- system.time(res1D <- run_sim(nt=1e4,rptfreq=10,seed=101,progress=TRUE,
                                    dt=0.2,hazard=TRUE))


hazdtvec <- 1/c(1,2,5,10,20)
hazdtList <- lapply(hazdtvec,
                    function(dt) {
    run_sim(nt=1e4,seed=101,progress=TRUE,
            dt=dt,hazard=TRUE)
    })

sizevec <- round(10^seq(log10(20),3,length.out=8))
sizeList <- lapply(sizevec,
                   function(N) {
    run_sim(mod_params=list(N=N),
            nt=1e5,seed=101,progress=TRUE,
            dt=20,hazard=TRUE)
    })

csizeList <- lapply(sizevec,
                   function(N) {
    run_sim(mod_params=list(N=N),
            nt=1e5,seed=101,progress=TRUE,
            discrete=FALSE)
    })

save_vars <- ls(pattern="((res|t)1[[:upper:]]|List$|vec$)")
save(list=save_vars,
     file="../simdata/nn11.rda")


