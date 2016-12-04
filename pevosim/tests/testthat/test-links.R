devtools::load_all("~/projects/nearlyneutralHP/pevosim")

ftest <- function(dt=1) {
    multcloglog(dt=dt)$linkinv(1)
}
ftest(1)
dtvec <- 10^seq(-1,1,length=101)
ff <- sapply(dtvec,ftest)
plot(dtvec,ff,type="l",log="x")
m <-     multcloglog(dt=2)
expect_equal(m$linkinv(m$linkfun(0.7)),0.7)
expect_equal(m$linkfun(m$linkinv(2)),2)
