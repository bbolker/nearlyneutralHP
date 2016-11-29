#include <Rcpp.h>
#include "pevosim.h"
#define MAXVAL 1.0e308


/*  was using:
// http://gallery.rcpp.org/articles/using-the-Rcpp-based-sample-implementation/
*/

using namespace Rcpp;

typedef std::vector<double> dvec;
typedef std::vector<int> ivec;

// [[Rcpp::export]]
int mySample(NumericVector rates, double sum_rates) {
    double val = runif(1)[0]*sum_rates;
    int w=0;
    while(val>=0) {
	val -= rates[w++];
    }
    return w;
}

// [[Rcpp::export]]
NumericVector get_ratesC(dvec betavec,
			 dvec gamma,
			 ivec Ivec,
			 int S,
			 bool &overflow) {

    int nstrains = betavec.size();

    // http://stackoverflow.com/questions/6399090/c-convert-vectorint-to-vectordouble
    dvec inf_rates(nstrains);
    dvec recover_rates(nstrains);

    overflow=false;
    for (int i=0; i<nstrains; i++) {
	inf_rates[i]=betavec[i]*Ivec[i]*S;
	recover_rates[i]= gamma[i]*Ivec[i];
	if (inf_rates[i]>MAXVAL || recover_rates[i]>MAXVAL) {
	    overflow=true;
	}
    }

    inf_rates.insert(inf_rates.end(), 
		     recover_rates.begin(),
		     recover_rates.end());
    return wrap(inf_rates);
}

void do_mutC(List state,
	    const String mut_var,
	    dvec orig_trait,
	    double mut_mean,
            double mut_sd) {

    // (for now ...) assume a *single* mutation;
    // assume exponential inverse-link
    dvec mutvec = as<dvec>(state[mut_var]);
    dvec ltraitvec = as<dvec >(state["ltraitvec"]);
    ivec Ivec = as<ivec >(state["Ivec"]);

    double mut_chg = rnorm(1,mut_mean,mut_sd)[0];
    double new_trait = orig_trait[0] + mut_chg;
    mutvec.push_back(exp(new_trait));
    ltraitvec.push_back(new_trait);
    Ivec.push_back(1);
    state[mut_var] = mutvec; // not sure if this is efficient ...
    state["ltraitvec"] = ltraitvec;
    state["Ivec"] = Ivec;
}

// direct mod version 
void do_mut2C(dvec &ltraitvec,
	     dvec &betavec,
	     dvec &gamma,
	     ivec &Ivec,
	     const String mut_var,
	     int strain,
	     double mut_mean,
	     double mut_sd,
             bool debug=false) {

    // (for now ...) assume a *single* mutation;
    // assume exponential inverse-link

    double new_trait;

    double orig_ltrait = ltraitvec[strain];
    double orig_beta = betavec[strain];
    double orig_gamma = gamma[strain];
    double mut_chg = rnorm(1,mut_mean,mut_sd)[0];
    double new_ltrait = orig_ltrait + mut_chg;

    if (debug) {
	int nstrain=betavec.size();
	Rprintf("%d/%d %lf %lf %lf %lf\n",
		strain,nstrain,
		orig_ltrait,mut_chg,new_trait,exp(new_trait));
    }
    
    if (mut_var=="beta") {
	betavec.push_back(exp(new_ltrait));
	gamma.push_back(orig_gamma);
    } else {
	new_trait = orig_gamma + mut_chg;
	gamma.push_back(exp(new_trait));
	betavec.push_back(exp(orig_beta));
    }
    ltraitvec.push_back(new_ltrait);
    Ivec.push_back(1);

}

// [[Rcpp::export]]
void do_extinctC(List state,
		const String mut_var,
		int extinct) {

    // make a *single* strain extinct

    dvec mutvec = as<dvec >(state[mut_var]);
    dvec ltraitvec = as<dvec >(state["ltraitvec"]);
    ivec Ivec = as<ivec >(state["Ivec"]);

    mutvec.erase(mutvec.begin() + (extinct-1));
    ltraitvec.erase(ltraitvec.begin() + (extinct-1));
    Ivec.erase(Ivec.begin() + (extinct-1));

    state[mut_var] = mutvec; // not sure if this is efficient ...
    state["ltraitvec"] = ltraitvec;
    state["Ivec"] = Ivec;
}

void do_extinct2(dvec &mutvec,
		 dvec &ltraitvec,
		 ivec &Ivec,
		 int extinct) {

    // make a *single* strain extinct

    // Rprintf("do_extinct2: %d\n",Ivec.size());
    mutvec.erase(mutvec.begin() + (extinct));
    ltraitvec.erase(ltraitvec.begin() + (extinct));
    Ivec.erase(Ivec.begin() + (extinct));
    // Rprintf("do_extinct2: %d\n",Ivec.size());

}

// [[Rcpp::export]]
void run_stepC(List state,double t_tot, double t_end,
	       List params, double maxrate=1000, bool debug=true) {

    if (debug) Rprintf("begin\n");
    dvec ltraitvec = as< dvec >(state["ltraitvec"]);
    ivec Ivec =      as< ivec >(state["Ivec"]);
    dvec betavec =   as< dvec >(state["beta"]);
    dvec gamma =     as< dvec >(state["gamma"]);
    int S = state["S"];
    if (debug) Rprintf("state loaded\n");

    // extract params: mu, mut_var, mut_link, mut_sd
    double mu = params["mu"];
    double mut_sd = params["mut_sd"];
    double mut_mean = params["mut_mean"];
    String mut_var = params["mut_var"];
    String mut_link = params["mut_link"];

    if (debug) Rprintf("params loaded\n");

    int nstrain, event, strain, w;
    double r_mut, tot_rates;
    NumericVector rates;
    bool overflow=false;

    while (t_tot<t_end) {
	nstrain = ltraitvec.size();	
	rates = get_ratesC(betavec,gamma,Ivec,S,overflow);
	if (overflow) break;
	if (debug) {
	    Rprintf("t=%lf n=%d S=%d I=%d minrate=%lf maxrate=%lf\n",t_tot,nstrain,
		    S,Ivec[0]);
	}
	tot_rates = sum(rates);
	t_tot += rexp(1,tot_rates)[0];
	if (debug) Rcout << "rates " << rates << std::endl;
	w = mySample(rates, tot_rates);
	
	event =  ((w-1) / nstrain) + 1;
	strain = ((w-1) % nstrain); // 0-index!
	if (debug) Rprintf("w: %d %d %d\n",w,event,strain);
	if (event==1) { // infection
	    if (debug) Rprintf("infection\n");
	    r_mut = runif(1)[0];

	    if (runif(1)[0]<mu) {
		if (debug) Rprintf("mutation\n");
		do_mut2C(ltraitvec,
			betavec,
			gamma,
			Ivec,
			mut_var,
			strain,
			mut_mean,
			mut_sd,
                        debug);
	    } else {
		Ivec[strain]++;
	    }
	    S--;
	} // infection
	if (event==2) { // recovery
	    if (debug) Rprintf("recovery\n");
	    Ivec[strain]--;
	    S++;
	    if (Ivec[strain]==0) {
		do_extinct2((mut_var=="beta" ? betavec : gamma),
			    ltraitvec,
			    Ivec,
			    strain);
		if (debug) Rprintf("remaining strains: %d\n",Ivec.size());
	    }
	    if (Ivec.size()==0) break;
	    // need to deal with complete extinction
	}
	// assertions:
	// stopifnot(length(S)==1)
	// stopifnot(sum(state$Ivec)+S == N)
    } // loop until end of time period
    if (debug) Rprintf("done\n");
    if (overflow) {
	Rf_error("Rate overflow");
    }
    state["ltraitvec"] = ltraitvec;
    state["Ivec"] = Ivec;
    state["beta"] = betavec;
    state["gamma"] = gamma;
    state["S"] = S;
}
