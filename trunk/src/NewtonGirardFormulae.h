/*
 * =====================================================================================
 *
 *       Filename:  NewtonGirardFormulae.h
 *
 *    Description:  Implementation of Newton-Girard formulae, which allows
 *									the conversion between power sum and elementary
 *									symmetric polynomial.
 *
 *        Version:  1.0
 *        Created:  8/26/2012 2:22:42 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu, hh1985@bu.edu
 *   Organization:  Boston University
 *
 * =====================================================================================
 */
#ifndef BRAIN_NEWTONGIRARDFORMULAE_H
#define BRAIN_NEWTONGIRARDFORMULAE_H

#include "Polynomials.h"

namespace msmath
{

	
	// Basic principles:
	// 1. The implementation of Newton-Girard formula should be independent of Vietes formula;
	// 2. The class is not supposed to store the values of power sums and ESPs.
  class NewtonGirardFormulae
  {
  private:
    // The number of variables, which is also equal to the number of roots in the polynomial equations.
    size_t _var_num;

		// Vector of power sums at different order number. _power_sum[0] arbitrarily set as 0 since it is not used for calculation.
    // PowerSumVec _power_sum;
    // EleSymPolyVec _ele_sym;

    // Updating power sum values given elementary symmetric polynomial.
    void updatePowerSum(PowerSumVec& ps_vec, const EleSymPolyVec& esp_vec);

    // Updating elementary symmetric polynomial given power sum values.
		void updateElementarySymmetricPoly(const PowerSumVec& ps_vec, EleSymPolyVec& esp_vec);
  
  public:
    // Constructor. 
		// NewtonGirardFormulae() {}

		// Preferred constructor.
		NewtonGirardFormulae(size_t var_num)
			: _var_num(var_num) {}

		inline void setOrder(size_t n)
		{
			_var_num = n;
		}
		inline double getOrder() const
		{
			return _var_num;
		}

		// This function will automatically fill the missing value between 
		// ps_vec and esp_vec.
		void updateParameters(PowerSumVec& ps_vec, EleSymPolyVec& esp_vec);
		inline void updateParameters(PolyParam& poly_param)
		{
			this->updateParameters(poly_param.second, poly_param.first);
		}

  };
}


#endif /* BRAIN_NEWTONGIRARDFORMULAE_H */