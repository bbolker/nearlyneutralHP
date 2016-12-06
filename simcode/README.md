# base runs

- `nn2`: basic runs for 10^6 steps for N=1000, 100, 33
- `nn3`: LHS. R0_init=4, 10^6 steps, N from 33 to 1000, variation in mu, mut_mean, mut_sd, gamma
- `nn4`: like `nn2` but with evolving gamma rather than beta
- `nn5`: like `nn3` but with evolving gamma rather than beta
- `nn6`: fix a few things from `nn5`
- `nn7`: ditto

# next steps

- `nn8`: Sobol sequence; continuous time/C++
- `nn12
