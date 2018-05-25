#include <R.h>
#include <Rmath.h>
 
#define NODEBUG

/****************************************************************************
 * Author: Yves Deville <deville.yves@alpestat.com>                         *
 *                                                                          *
 * The 'Scale' function computes the values of the scaling at the design    *
 * points                                                                   *
 *                                                                          *
 *                                                                          *
 ****************************************************************************/

void Scale(int *n,               // length of x
	   int *nKnots,          // length of knots
	   int *icuts,           // array of cuts
	   double *x,            // design points
	   double *knots,        // knots 
	   double *eta,          // fun values at knots
	   double *scale) {      // function values   
  
  int 
    i, ell;

  double 
    S,
    knotL, knotR, 
    dKnot,
    etaL, etaR,
    deltaL, deltaR, diffL;

  /*========================================================
   * the knots define (nKnots - 1) intervals numeroted 
   * from 0 to nKnots -1
   *=======================================================*/
  
  S = 0.0;

  for (ell = 0; ell < *nKnots -1; ell++) {
 
    // interval ell
    knotR = knots[ell + 1];
    knotL = knots[ell];
    dKnot =  knotR - knotL;

    etaL = eta[ell];
    etaR = eta[ell+1];
    
    for (i = icuts[ell]; i < icuts[ell+1]; i++) {
      
      deltaL = x[i] - knotL;
      deltaR = knotR - x[i];
      diffL = deltaL / dKnot;
      
      scale[ i ] = S + 0.5 * diffL * (etaL * (dKnot + deltaR) + etaR*deltaL);
      
    }
       
    S += 0.5 * ( etaL + etaR ) * dKnot; 

  }
	 
}

/****************************************************************************
 * Compute scaling 'gradients', i.e. derivative w.r.t. the parameters 'eta' *
 *                                                                          *
 *                                                                          *
 ****************************************************************************/

void gradScale(int *n,               // length of x
	       int *nKnots,          // length of knots
	       int *icuts,           // array of cuts
	       double *x,            // design points
	       double *knots,        // knots 
	       double *grad) {       // gradient array   
  
  int 
    i, ell;

  double 
    knotL, knotR, knot, 
    dKnotL, dKnotR, dKnot,
    dif;

  /*========================================================
   * the knots define (nKnots - 1) intervals numeroted 
   * from 0 to nKnots -1
   *=======================================================*/
 
  for (ell = 0; ell < *nKnots; ell++) {

    // interval ell-1, if it exists
    if (ell > 0) {

      knotL = knots[ell - 1];
      dKnotL = knots[ell] - knotL;

      for (i = icuts[ell-1]; i < icuts[ell]; i++) {
	dif = x[i] - knotL;
	grad[i + ell * *n ] = 0.5 * dif * dif / dKnotL;
      }

    } else {
      dKnotL = 0.0;
    }

    if (ell < *nKnots -1) {
      // interval ell
      knotR = knots[ell + 1];
      dKnotR = knotR - knots[ell];
      dKnot =  dKnotL + dKnotR;
      
      for (i = icuts[ell]; i < icuts[ell+1]; i++) {
	dif = knotR - x[i];
	grad[i + ell * *n ] = 0.5* (dKnot - dif * dif / dKnotR);
      }
      
      //intervals >= ell + 1, if any 
      if (ell < *nKnots -2) {
	
	for (i = icuts[ell+1]; i < *n; i++) {
	  grad[i + ell * *n ] = 0.5 * dKnot;
	}
	
      } 
    }

  }
	 
}
