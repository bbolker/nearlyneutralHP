// http://gallery.rcpp.org/articles/using-the-Rcpp-based-sample-implementation/
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_rates(std::vector<double> betavec,
			std::vector<double> gamma,
			std::vector<int> Ivec,
			int S) {

    int nstrains = betavec.size();
    int nrates = 2*nstrains;
    int i;
    NumericVector rates(nrates);
    // http://stackoverflow.com/questions/6399090/c-convert-vectorint-to-vectordouble
    std::vector<int> dIvec(Ivec.begin(), Ivec.end());
    std::vector<double> inf_rates(nstrains);
    std::vector<double> recover_rates(nstrains);

    for (i=0; i<nstrains; i++) {
	inf_rates[i]=betavec[i]*dIvec[i]*S;
	recover_rates[i]= gamma[i]*dIvec[i];
    }
    inf_rates.insert(inf_rates.end(), 
		     recover_rates.begin(),
		     recover_rates.end());
    return wrap(inf_rates);
}

// void test_list(List state, String ind) {
//    state[ind] = 1.0;
// }
    

List mk_list(double x, double y) {
    List L(2);
    std::vector<double> v1(1), v2(1);
    v1[0] = x;
    v2[0] = y;
    L[1] = v1;
    L[2] = v2;
    return L;
}

// [[Rcpp::export]]

void do_mut(List state,
	    const String mut_var,
	    std::vector<double> orig_trait,
	    double mut_mean,
	    double mut_sd) {

    // (for now ...) assume a *single* mutation;
    // assume exponential inverse-link
    std::vector<double> mutvec = as<std::vector<double> >(state[mut_var]);
    std::vector<double> ltraitvec = as<std::vector<double> >(state["ltraitvec"]);
    std::vector<int> Ivec = as<std::vector<int> >(state["Ivec"]);

    double mut_chg = rnorm(1,mut_mean,mut_sd)[0];
    double new_trait = orig_trait[0] + mut_chg;
    mutvec.push_back(exp(new_trait));
    ltraitvec.push_back(new_trait);
    Ivec.push_back(1);
    state[mut_var] = mutvec; // not sure if this is efficient ...
    state["ltraitvec"] = ltraitvec;
    state["Ivec"] = Ivec;
}

// [[Rcpp::export]]

void do_extinct(List state,
		const String mut_var,
		int extinct) {

    // make a *single* strain extinct

    std::vector<double> mutvec = as<std::vector<double> >(state[mut_var]);
    std::vector<double> ltraitvec = as<std::vector<double> >(state["ltraitvec"]);
    std::vector<int> Ivec = as<std::vector<int> >(state["Ivec"]);

    mutvec.erase(mutvec.begin() + (extinct-1));
    ltraitvec.erase(ltraitvec.begin() + (extinct-1));
    Ivec.erase(Ivec.begin() + (extinct-1));

    state[mut_var] = mutvec; // not sure if this is efficient ...
    state["ltraitvec"] = ltraitvec;
    state["Ivec"] = Ivec;
}

// [[Rcpp::export]]

void run_stepC(List state,double t_tot, double t_end,
	       List params, bool debug=true) {

    if (debug) Rprintf("begin\n");
    std::vector<double> ltraitvec = 
	as< std::vector<double> >(state["ltraitvec"]);
    std::vector<int> Ivec = 
	as< std::vector<int> >(state["Ivec"]);
    std::vector<double> betavec = 
	as< std::vector<double> >(state["beta"]);
    std::vector<double> gamma = 
	as< std::vector<double> >(state["gamma"]);
    int S = state["S"];
    if (debug) Rprintf("state loaded\n");

    // extract params: mu, mut_var, mut_link, mut_sd
    double mu = params["mu"];
    double mut_sd = params["mut_sd"];
    double mut_mean = params["mut_mean"];
    String mut_var = params["mut_var"];
    String mut_link = params["mut_link"];

    if (debug) Rprintf("params loaded\n");

    int nstrain = ltraitvec.size();	
    int nrates = 2*nstrain;
    int event, strain;
    double r_mut;
    NumericVector rates(nrates);
    IntegerVector wseq = seq(1,nrates);
    int w;

    while (t_tot<t_end) {
	nstrain = ltraitvec.size();	
	rates = get_rates(betavec,gamma,Ivec,S);
	if (debug) {
	    Rprintf("t=%lf %d %lf %lf\n",t_tot,nstrain,
			   min(rates),max(rates));
	}
	t_tot += rexp(1,sum(rates))[0];
	
	w = RcppArmadillo::sample(wseq, nrates, false, rates)[0];
	// w = sample(length(rates),size=1,prob=rates);

	event =  ((w-1) / nstrain) + 1;
	strain = ((w-1) % nstrain) + 1;
	if (event==2) { // recovery
	    Ivec[strain]--;
	    S++;
	    if (Ivec[strain]==0) {
		do_extinct(state,mut_var,strain);
	    }
	}
	if (event==1) { // infection
	    if (debug) Rprintf("infection\n");
	    r_mut = runif(1)[0];
	    std::vector<double> orig_ltrait(1);
	    orig_ltrait[0] = ltraitvec[strain];

	    if (runif(1)[0]<mu) {
		do_mut(state,
		       mut_var,
		       // mut_link,
		       orig_ltrait,
		       mut_mean,
		       mut_sd);
	    } else {
		Ivec[strain]++;
	    }
	    S--;
	} // infection
	// assertions:
	// stopifnot(length(S)==1)
	// stopifnot(sum(state$Ivec)+S == N)
    } // loop until end of time period
    if (debug) Rprintf("done\n");
    state["S"] = S;
}
