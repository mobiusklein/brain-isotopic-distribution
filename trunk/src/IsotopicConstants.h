/********************************************************************
	created:	2012/08/24
	created:	24:8:2012   10:23
	filename: 	IsotopicConstants.h
	file path:	GAGPL\SPECTRUM
	file base:	IsotopicConstants
	file ext:	h
	author:		Han Hu
	
	purpose:	It is used to store the temporary element constants values								used to calculate theoretical isotopic distribution. Only one 
						copy is maintained.
*********************************************************************/
#ifndef BRAIN_ISOTOPICCONSTANTS_H 
#define BRAIN_ISOTOPICCONSTANTS_H

#include "Element.h"
#include "NewtonGirardFormulae.h"
#include "VietesFormulae.h"
#include "Polynomials.h"
#include "PeriodicTable.h"
#include <boost/noncopyable.hpp>
#include <map>

namespace brain 
{
	using namespace msmath;

	// The element constants and modified elemental constants (mass added).
	//typedef std::pair<PolyParam, PolyParam> PhiConstants;
	struct PhiConstants
	{
		size_t nth_order;
		PolyParam ele_param;
		PolyParam ele_mass_param;
	};

	// The element symbol and corresponding constants
	typedef std::map<std::string, PhiConstants> ElementConstants;

	class IsotopicConstants : private boost::noncopyable
	{
	private:
		// The map stores the basic parameters for each element.
		ElementConstants _ele_const;

		// The map stores the high order for each element. The information can be freed if needed.
		//ElementConstants _ele_high_order;

		// Number of variants arbitrarily defined.
		size_t _var_num;

		PeriodicTable& ptable;

	private:
		IsotopicConstants(size_t order = 0)
			: _var_num(order), ptable(PeriodicTable::Instance()) {}
		
		// If has_mass == true, it means the parameters for calculating 
		// centered mass. This function is only element-related, and it has 
		// nothing to do with _var_num.
		PolyParam updateCoefsOfElement(const Element& e, bool has_mass = false);

		// Iterate over all elements, update their constants to the order of _var_num.
		void update();

	public:
		// Singleton. Set the initial value as 0.
		static IsotopicConstants& Instance()
		{
			static IsotopicConstants iso_const;
			return iso_const;
		}
		
		// Updating the _var_num to specified order, if the specified order is smaller than _var_num, there is no need to update.
		void updateOrder(size_t order);

		// Cleaning the constants.
		inline void clear()
		{
			_ele_const.clear();
		}

		// Getting the constants for calculating probabilities.
		inline const ElementConstants& getProbConstants()
		{
			return _ele_const;
		}

		// (r^-1)^l for specified element.
		double getNthElementPowerSum(const std::string& symbol, size_t order);

		// (r^-1)^l for mass-modified element formula
		double getNthModifiedElementPowerSum(const std::string& symbol, size_t order);

		// Adding the most basic information into the _ele_const map.
		// This function has nothing to do with the _var_num.
		// The basic assumption is that once one element has been stored, the minimum of order is specified.
		void addElement(const std::string& symbol);
	};
}

#endif /* BRAIN_ISOTOPICCONSTANTS_H */