# base runs

- `nn2`: basic runs for 10^6 steps for N=1000, 100, 33
- `nn3`: LHS. R0_init=4, 10^6 steps, N from 33 to 1000, variation in mu, mut_mean, mut_sd, gamma
- `nn4`: like `nn2` but with evolving gamma rather than beta
- `nn5`: like `nn3` but with evolving gamma rather than beta
- `nn6`: fix a few things from `nn5`
- `nn7`: ditto


- `nn8`: Sobol sequence; continuous time/C++

- `nn11`: basic runs for CMS presentation
   - 1A-D: all 10^4 steps, discrete, continuous, hazard with dt=1 and 0.2
   - hazdt: hazard with dt=1/(c(1,2,5,10,20))
   - size: pop size from 20 to 1000. dt=20 -- should have been 1/20!
   - csize: ditto, continuous-time
