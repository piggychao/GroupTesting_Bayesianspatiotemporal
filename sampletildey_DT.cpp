#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// sampletildey_DT - Implements posterior sampling of y_tilde under the Dorfman testing (DT) protocol
//
// INPUT:
// pooli     - 1×n vector of pool indices.
// pind     - 1×n vector of individual probabilities.
// z        - 1×J vector of pool responses.
// yinit    - 1×n vector of initial y_tilde values.
// yretest  - 1×n vector of individual retesting results.  
// Se       - 1×2 vector of sensitivities for pool and individual tests.  
// Sp       - 1×2 vector of specificities for pool and individual tests.  
//
// OUTPUT:
// yinit    - 1×n vector of posterior samples for y_tilde.
//
// [[Rcpp::export]]
NumericVector sampletildey_DT(NumericVector pooli, 
                            NumericVector pind, NumericVector z, 
                            NumericVector yinit, NumericVector yretest, 
                            NumericVector Se, NumericVector Sp){
  double nsample = yinit.length();
  double npool = z.length();
  int size = ceil(nsample / npool);
  int endid;
  double Se1, Se2, Sp1, Sp2;
  Se1 = Se[0]; //pool sensitivity
  Se2 = Se[1]; // individual testing sensitivity
  Sp1 = Sp[0]; // pool specificity
  Sp2 = Sp[1]; // individual testing specificity
  
  for(int i=0; i<nsample; ++i){
    double pind0, p1, p0, pnew;
    int id = pooli[i]; // get pool index for individual i
    int startid = id*size;
    endid = id*size+(size-1);
    if (endid>(nsample-1)) {endid=nsample-1;}
    IntegerVector poolind = seq(startid, endid); // use pool-individual ids to find subset pool
    pind0 = pind[i]; // probability of individual i
    double z0 = z[id]; // get the corresponding pool response
    NumericVector spool = yinit[poolind]; // find corresponding subset pool
    int yind0 = i - size*id; // get individual index in subset pool
    spool.erase(spool.begin() + yind0); // removing y_i from subset pool
    double spools = std::accumulate(spool.begin(), spool.end(), 0); // sub-pool result after removing i
    double z1;
    if (spools > 0) {z1=1;} else {z1=0;} // pool response for new sub-pool
    // Under the DT protocol
    if (z0 > 0){
      double retest0 = yretest[i];
      p1 = pind0 * pow(Se1, z0) * pow(1-Se1, (1-z0)) *
        pow(Se2, retest0) * pow(1-Se2, 1-retest0);
      p0 = (1-pind0) * pow( pow(Se1, z0)*pow(1-Se1, 1-z0), z1) *
        pow( pow(1-Sp1, z0)*pow(Sp1,1-z0), 1-z1) *
        pow(1-Sp2, retest0) * pow(Sp2, 1-retest0);
    } else {
      p1 = pind0 * pow(Se1, z0) * pow(1-Se1, (1-z0));
      p0 = (1-pind0) * pow( pow(Se1, z0)*pow(1-Se1, 1-z0), z1) *
      pow( pow(1-Sp1, z0)*pow(Sp1, 1-z0), 1-z1); 
    }
    pnew = p1/(p1+p0);
    yinit[i] = R::rbinom(1, pnew);
  }
  return(yinit);
}


