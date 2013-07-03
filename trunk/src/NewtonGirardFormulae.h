/*
 * =====================================================================================
 *
 *       Filename:  NewtonGirardFormulae.h
 *
 *    Description:  
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

#include <vector>
#include "VietesFormulae.h"

namespace msmath
{
	typedef std::vector<double> PowerSumVec;
	

  class NewtonGirardFormulae
  {
  private:
    // The highest order of the polynomial is determined by _coef.size()-1.
    //const PolyCoef& _coef; 
    size_t _var_num;
    PowerSumVec _power_sum;
    EleSymPolyVec _ele_sym;

    // Update power sum values given elementary symmetric polynomial.
		// If empty, just start from the scratch.
		// Otherwise, build from the current last element.
    void updatePowerSum();

    // Similarly, update elementary symmetric polynomial given power sum values
		void updateElementarySymmetricPoly();
  
  public:
    // Constructor. 
		NewtonGirardFormulae() {}
		NewtonGirardFormulae(size_t var_num)
			: _var_num(var_num) {}

		inline void setOrder(size_t n)
		{
			_var_num = n;
		}
		inline double getOrder() 
		{
			return _var_num;
		}

		void setElementarySymmetricPoly(const EleSymPolyVec& ele_vec);

		void setPowerSum(const PowerSumVec& ps_vec);

		inline PowerSumVec getPowerSumVec()
		{
			return _power_sum;
		}
		inline EleSymPolyVec getElementarySymmetricPolyVec()
		{
			return _ele_sym;
		}
		inline double getPowerSum(const size_t index)
		{
			return _power_sum.at(index);
		}
		inline double getElementarySymmetricPoly(const size_t index)
		{
			return _ele_sym.at(index);
		}

		void updateToHigherOrder(const size_t k_th_order);

		// The constriction of the definition of e is implemented by this function.
		void updateUsingVietesFormulae(PolyCoef& coef, const size_t k_th_order);
		void updateUsingVietesFormulae(PolyCoef& coef);


  };
}


#endif /* BRAIN_NEWTONGIRARDFORMULAE_H */