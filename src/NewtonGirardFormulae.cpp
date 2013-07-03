#include "NewtonGirardFormulae.h"

namespace msmath
{
  
	void NewtonGirardFormulae::updateElementarySymmetricPoly()
	{
		if(_ele_sym.empty()) {
			// e(0) = 1
			_ele_sym.push_back(1);
		}
		// Matching elementary symmetric polynomial values using power sum values.
		if(_ele_sym.size() >= _power_sum.size()) {
			return;
		}


		for(size_t i=_ele_sym.size(); i<_power_sum.size(); i++)
		{
			double temp_ele = 0.0;
			
			for(size_t j=1; j<=i; j++)
			{
				int sign = (j%2 == 1 ? 1 : -1);
				temp_ele += sign * _power_sum.at(j) * _ele_sym.at(i-j); 
			}
			
			temp_ele /= i;
			
			_ele_sym.push_back(temp_ele);
		}

	}
  
	void NewtonGirardFormulae::updatePowerSum()
	{
	
		if(_power_sum.empty()) {
			// p(0) = 0. The value might be changed in the middle of calculation.
			_power_sum.push_back(0);
		}

		if(_power_sum.size() >= _ele_sym.size()) {
			return;
		}

		for(size_t i=_power_sum.size(); i<_ele_sym.size(); i++)
		{
			double temp_ps = 0.0;

			for(size_t j=1; j<=i; j++)
			{
				// Temporarily change the value of p(0).
				_power_sum[0] = j;
				int sign = (j%2 == 1 ? 1 : -1);
				temp_ps += sign * _ele_sym.at(j) * _power_sum.at(i-j);
			}

			_power_sum.push_back(temp_ps);
		}
		_power_sum[0] = 0;
	}

	
	void NewtonGirardFormulae::updateUsingVietesFormulae(PolyCoef& coef, const size_t k_th_order)
	{
		// update _var_num.
		_var_num = coef.size()-1;

		_ele_sym.clear();
		
		// e(0) = 1.
		_ele_sym.push_back(1);

		VietesFormulae formulae(coef);
		
		EleSymPolyVec temp_ele_vec = formulae.getElementarySymmetricFunctionFromCoef();

		_ele_sym.insert(_ele_sym.end(), temp_ele_vec.begin(),temp_ele_vec.end());

		//updateToHigherOrder(k_th_order);
		updatePowerSum();

	}

	void NewtonGirardFormulae::updateUsingVietesFormulae(PolyCoef& coef )
	{
		const size_t k = coef.size()-1;
		this->updateUsingVietesFormulae(coef, k);
	}

	void NewtonGirardFormulae::updateToHigherOrder(const size_t k_th_order )
	{
		if(k_th_order > _var_num) {
			//_var_num = k_th_order;
			// Update corresponding power sum values.
			for(size_t i=_var_num+1; i<=k_th_order; i++)
			{
				// Update elementary symmetric functions.
				_ele_sym.push_back(0);
			}
			// Update corresponding power sum values.
			this->updatePowerSum();
		} else {
			// Do nothing:)
		} 
	}

	void NewtonGirardFormulae::setElementarySymmetricPoly(const EleSymPolyVec& ele_vec)
	{
		_ele_sym = ele_vec;
		_power_sum.clear();
		this->updatePowerSum();
	}

	void NewtonGirardFormulae::setPowerSum(const PowerSumVec& ps_vec)
	{
		_power_sum = ps_vec;
		// The update function will not clean the values. 
		_ele_sym.clear();
		this->updateElementarySymmetricPoly();
	}



}