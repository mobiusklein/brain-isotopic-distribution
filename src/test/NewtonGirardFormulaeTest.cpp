/*
 * =====================================================================================
 *
 *       Filename:  NewtonGirardFormulaeTest.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  8/26/2012 5:12:38 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu, hh1985@bu.edu 
 *   Organization:  Boston University.
 *
 * =====================================================================================
 */
#include "NewtonGirardFormulae.h"
#include <iostream>

int main() {
  // Create a vector of coefficient for n-order polynomial equation.

	using namespace msmath;
	using namespace std;

	double mycoef[] = {3, 2, 0, 5, 1};
  PolyCoef poly_eq(mycoef, mycoef + sizeof(mycoef) / sizeof(double));

  NewtonGirardFormulae formulae;
	formulae.updateUsingVietesFormulae(poly_eq);
	formulae.updateToHigherOrder(7);

  for(int i=0; i<7; i++) {
    cout << "ID\tPowersum\tElementary symmetric polynomial" << endl;
    cout << i << "\t" << formulae.getPowerSum(i) << "\t" << formulae.getElementarySymmetricPoly(i) << endl;
  }

	cin.get();
}