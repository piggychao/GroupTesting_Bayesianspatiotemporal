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
// sampletildey_AT - Implements posterior sampling of y_tilde under the array testing (AT) protocol
//
// INPUT:
// pooli    - n×2 matrix of pool indices (col1 - horizontal pool indices; col2 - vertical pool indices).
// width    - Array width
// vpooli   - (Array width × number of vertical pools) matrix of individual indices in vertical pools.
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
//
// [[Rcpp::export]]
NumericVector sampletildey_AT(NumericMatrix pooli, IntegerMatrix vpooli, double width,
                            NumericVector pind,
                   NumericVector z, 
                   NumericVector yinit, NumericVector yretest,
                   NumericVector Se, NumericVector Sp) {
 
  int startid, endid, poolid0, yind0;
  double spools, z1, p1, p0, pnew, pind0, poolp1, poolp0; // sum of subpool
  double nsample = yinit.length();
  //int npool = z.length();
  int size = width;
  int nspool = pooli.ncol();
  double nhpool = ceil(nsample/width);
  IntegerVector poolind;
  NumericVector id, spool, z0v(nspool), z1v(nspool); // horizontal and vertical pools for individual
  LogicalVector yindlogic, pindlogic;
  double Se1, Se2, Sp1, Sp2;
  Se1 = Se[0]; //pool sensitivity
  Se2 = Se[1]; // individual testing sensitivity
  Sp1 = Sp[0]; // pool specificity
  Sp2 = Sp[1]; // individual testing specificity
  
  for(int i=0; i<nsample; ++i){
   id = pooli(i, _); // get pool index for individual i
   pind0 = pind[i]; // probability of individual i
   poolp1 = pind0;
   poolp0 = 1-pind0;
    
   for(int j=0; j<nspool; ++j){
     poolid0 = id[j];
     z0v[j] = z[poolid0];

     if (poolid0 < nhpool){
        startid = poolid0*size;
        endid = poolid0*size+(size-1);
        if (endid >= nsample) {endid = nsample-1;} 
        // use pool-individual ids to find subset pool
        poolind = seq(startid, endid);
        spool = yinit[poolind];
        yind0 = i - size*poolid0; // get individual index in subset pool
        spool.erase(spool.begin() + yind0); // removing y_i from subset pool
        spools = std::accumulate(spool.begin(), spool.end(), 0);
        if (spools > 0) {z1=1;} else {z1=0;}
      } else {
        poolind = vpooli(poolid0 - nhpool, _);
        pindlogic = !is_na(poolind);
        poolind = poolind[pindlogic];
        spool = yinit[poolind];
        yindlogic = !(poolind == i); // get individual logic index in subset pool
        spool = spool[yindlogic]; // removing y_i from subset pool
        spools = std::accumulate(spool.begin(), spool.end(), 0);
        if (spools > 0) {z1=1;} else {z1=0;}
      }
      z1v[j] = z1;
      
      // Under the AT protocol
      poolp1 = poolp1 * pow(Se1, z0v[j]) * pow(1-Se1, 1-z0v[j]);
      poolp0 = poolp0 * pow(pow(Se1, z0v[j])*pow(1-Se1, 1-z0v[j]), z1) * 
        pow( pow(1-Sp1, z0v[j])*pow(Sp1,1-z0v[j]), 1-z1);
   }
  
  NumericVector poolres = z[id];
  if (std::accumulate(poolres.begin(), 
                      poolres.end(), 0) == nspool) 
  {double retest0 = yretest[i];
   p1 = poolp1 * pow(Se2, retest0) * pow(1-Se2, 1-retest0);
   p0 = poolp0 * pow(1-Sp2, retest0) * pow(Sp2, 1-retest0);}
  else {
   p1 = poolp1;
   p0 = poolp0;
  }
  pnew = p1/(p1+p0);
 // Rprintf("the value of pnew : %f \n", pnew);
  yinit[i] = R::rbinom(1, pnew);
  }
  return yinit;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

