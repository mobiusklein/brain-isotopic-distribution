/********************************************************************
	created:	2012/08/27
	created:	27:8:2012   14:30
	filename: 	VietesFormulae.cpp
	file path:	GAGPL\MATH
	file base:	VietesFormulae
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include "VietesFormulae.h"

namespace msmath
{
	EleSymPolyVec VietesFormulae::getElementarySymmetricFunctionFromCoef()
	{
		EleSymPolyVec ele_vec;
		
		// Highest order of the polynomial.
		const size_t max_n = _coef.size()-1;

		// e0 = 1 based on Newton-Girard formulae.
		// ele_vec.push_back(1);

		// Start from a(n).
		//double a_n = coefficients.back();
		for(size_t i = 1; i<=max_n; i++)
		{
			int sign = (i%2 == 0 ? 1 : -1);
			double ele_value = sign * _coef.at(max_n - i) / _coef.back(); 
			ele_vec.push_back(ele_value);
		}

		return ele_vec;
	}

}