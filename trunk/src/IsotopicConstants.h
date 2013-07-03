/********************************************************************
	created:	2012/08/24
	created:	24:8:2012   10:23
	filename: 	IsotopicConstants.h
	file path:	GAGPL\SPECTRUM
	file base:	IsotopicConstants
	file ext:	h
	author:		Han Hu
	
	purpose:	
*********************************************************************/
#ifndef BRAIN_ISOTOPICCONSTANTS_H 
#define BRAIN_ISOTOPICCONSTANTS_H

#include "Element.h"
#include "NewtonGirardFormulae.h"
#include <boost/noncopyable.hpp>
#include <map>

namespace brain 
{
	// Phi factor.
	using namespace msmath;
	// The element symbol and corresponding constants 
	typedef std::pair<NewtonGirardFormulae, NewtonGirardFormulae> PhiConstants;
	typedef std::map<std::string, PhiConstants> ElementConstants;

	class IsotopicConstants : private boost::noncopyable
	{
	private:
		// phi_vec[0] should be set as 0.
		ElementConstants _ele_const;

		size_t _var_num;

	private:
		IsotopicConstants()
			: _var_num(5) {}
		
		NewtonGirardFormulae updateCoefsOfElement(const Element& e, bool has_mass = false);

		// Iterate over all elements, update them based on _var_num.
		void update();

	public:
		// Singleton.
		static IsotopicConstants& instance();
		
		inline void setOrder(size_t order)
		{
			_var_num = order;
		}
		void updateOrder(size_t order);

		inline const ElementConstants& getProbConstants()
		{
			return _ele_const;
		}
		// (r^-1)^l for specified element. 
		inline double getNthPhiFactor(const std::string& symbol, size_t order)
		{
			return _ele_const[symbol].first.getPowerSum(order);
		}

		inline double getNthModifiedPhiFactor(const std::string& symbol, size_t order)
		{
			return _ele_const[symbol].second.getPowerSum(order);
		}

		void addElement(const std::string& symbol);
	};
}

#endif /* BRAIN_ISOTOPICCONSTANTS_H */