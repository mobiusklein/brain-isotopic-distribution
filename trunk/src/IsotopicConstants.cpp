#include "IsotopicConstants.h"
#include "PeriodicTable.h"
#include "NewtonGirardFormulae.h"
#include <numeric>

namespace brain 
{
	IsotopicConstants& IsotopicConstants::instance()
	{
		static IsotopicConstants iso_const;
		return iso_const;
	}

	NewtonGirardFormulae IsotopicConstants::updateCoefsOfElement( const Element& e, bool has_mass /*= false*/ )
	{
		// Get the coefficients of given element
		IsotopeSetSequential::reverse_iterator it = e.isotopes.rbegin();
		// Get the highest isotope neutron shift.
		size_t max_n = it->neutron_shift;

		std::vector<double> temp_vec;

		for(; it != e.isotopes.rend(); it++)
		{
			size_t current_order = max_n - it->neutron_shift;

			double coef = (has_mass == true ? it->mass : 1.0);

			// Empty shift.
			if(current_order > temp_vec.size()) {

				for(size_t i=temp_vec.size(); i<current_order; i++)
				{
					temp_vec.push_back(0.0);					
				}
				temp_vec.push_back(it->abundance * coef);

			} else if(current_order == temp_vec.size()){

				temp_vec.push_back(it->abundance * coef);		
			} else {

				throw std::runtime_error("The list of neutron shift is not ordered.");
			}

		}

		NewtonGirardFormulae newton;

		newton.updateUsingVietesFormulae(temp_vec);
		newton.updateToHigherOrder(_var_num);		

		return newton;
	}

	void IsotopicConstants::updateOrder( size_t order )
	{
		_var_num = order;
		update();
	}

	void IsotopicConstants::update()
	{
		ElementConstants::iterator it = _ele_const.begin();
		for(; it != _ele_const.end(); it++)
		{
			it->second.first.updateToHigherOrder(_var_num);
			it->second.second.updateToHigherOrder(_var_num);
		}
	}

	void IsotopicConstants::addElement(const std::string& symbol )
	{
		// If element has been added, ignore the rest part.
		ElementConstants::iterator it = _ele_const.find(symbol);
		if(it!=_ele_const.end())
			return;

		PeriodicTable& ptable = PeriodicTable::Instance();
		//ptable.load();
		Element ele = ptable.getElementBySymbol(symbol);

		PhiConstants phi_const;
		phi_const.first = this->updateCoefsOfElement(ele);
		phi_const.second = this->updateCoefsOfElement(ele, true);

		_ele_const.insert(std::make_pair(symbol, phi_const));
	}
}
