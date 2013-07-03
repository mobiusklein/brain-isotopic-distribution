/********************************************************************
	created:	2012/08/27
	created:	27:8:2012   8:40
	filename: 	E:\Sync\Dropbox\Codes\glycan_pipeline\GAG\src\GAGPL\SPECTRUM\IsotopicDistribution.cpp
	file path:	E:\Sync\Dropbox\Codes\glycan_pipeline\GAG\src\GAGPL\SPECTRUM
	file base:	IsotopicDistribution
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/
#include "IsotopicDistribution.h"
#include "Param.h"
#include <boost/make_shared.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/lexical_cast.hpp>
#include <algorithm>
#include <sstream>
#include <fstream>


namespace brain
{
	void IsotopicDistribution::createMonoPeak()
	{
		// Get the monoisotopic isotope for each element of the composition.
		monopeak.mz = _compo.getMass();

		const std::map<std::string, int> ele_count = _compo.get();

		std::map<std::string, int>::const_iterator it = ele_count.begin();

		for(; it != ele_count.end(); it++)
		{
			// For each element, get the value of prob^coef.
			PeriodicTable& ptable = PeriodicTable::Instance();
			Isotope iso = ptable.getIsotopeByRelativeShift(it->first);
			monopeak.intensity += it->second * log(iso.abundance);
		}
		monopeak.intensity = exp(monopeak.intensity);

	}
	
	double IsotopicDistribution::calculatePhiValue(size_t order)
	{
		double phi_value = 0.0;
		// For each element, get the phi factor.
		const std::map<std::string, int> ele_count = _compo.get();

		std::map<std::string, int>::const_iterator it = ele_count.begin();

		for(; it != ele_count.end(); it++)
		{
			phi_value += _iso_const.getNthPhiFactor(it->first, order) * it->second;
		}

		return phi_value;
	}

	PowerSumVec IsotopicDistribution::calculatePhiValueVec()
	{
		PowerSumVec phi_value_vec;
		// phi(0) = 0.
		phi_value_vec.push_back(0);

		for(size_t i=1; i<=_var_num; i++)
		{
			phi_value_vec.push_back(this->calculatePhiValue(i));
		}
			
		return phi_value_vec;
	}

	EleSymPolyVec IsotopicDistribution::calculateProbabilityVec()
	{
		PowerSumVec phi_vec = this->calculatePhiValueVec();
		
		NewtonGirardFormulae newton;
		newton.setPowerSum(phi_vec);

		EleSymPolyVec prob_vec = newton.getElementarySymmetricPolyVec();
	
		// q(j) = q(0) * e(j) * (-1)^j.
		for(size_t i=0; i<prob_vec.size(); i++)
		{
			int sign = (i%2 == 0 ? 1 : -1);
			prob_vec[i] = prob_vec[i] * monopeak.intensity * sign;
		}

		return prob_vec;
	}

	double IsotopicDistribution::calculateModifiedPhiValue(const std::string& symbol, size_t order)
	{
		double phi_value = 0.0;
		// For each element, get the phi factor.
		const std::map<std::string, int> ele_count = _compo.get();

		std::map<std::string, int>::const_iterator it = ele_count.begin();

		for(; it != ele_count.end(); it++)
		{
			int coef = (it->first == symbol ? it->second - 1 : it->second);
			phi_value += _iso_const.getNthPhiFactor(it->first, order) * coef;
		}

		phi_value += _iso_const.getNthModifiedPhiFactor(symbol, order); 

		return phi_value;
	}

	PowerSumVec IsotopicDistribution::calculateModifiedPhiValueVec(const std::string& symbol)
	{
		PowerSumVec phi_value_vec;
		// phi(0) = 1.
		phi_value_vec.push_back(0);

		for(size_t i=1; i<=_var_num; i++)
		{
			phi_value_vec.push_back(this->calculateModifiedPhiValue(symbol, i));
		}

		return phi_value_vec;
	}

	std::vector<double> IsotopicDistribution::calculateCenterMassVec(const EleSymPolyVec& prob_vec)
	{
		std::vector<double> mass_vec;

		const std::map<std::string, int> ele_count = _compo.get();

		PeriodicTable& ptable = PeriodicTable::Instance();

		NewtonGirardFormulae newton;
		
		for(size_t i=0; i<=_var_num; i++)
		{
			std::map<std::string, int>::const_iterator it = ele_count.begin();
			
			int sign = (i%2 == 0 ? 1 : -1);

			double center_mass = 0.0;

			for(; it != ele_count.end(); it++)
			{				
				const Isotope& mono = ptable.getIsotopeByRelativeShift(it->first);
				
				newton.setPowerSum(this->calculateModifiedPhiValueVec(it->first));
				
				// For each element, multiply the number of atoms by the modified phi.
				center_mass += it->second * sign * newton.getElementarySymmetricPoly(i) * monopeak.intensity * mono.mass;

			}
			// m(j) = sum(m(jk) * p(jk))/sum(p(jk))
			mass_vec.push_back(center_mass/prob_vec.at(i));
		}

		return mass_vec;
	}

	void IsotopicDistribution::setOrder(int order)
	{
		if(order == -1)
			_var_num = _compo.getMaxNumVariants();
		else
			_var_num = (size_t)order;
		_iso_const.updateOrder(_var_num);
	}

	AggregatedIsotopicVariants IsotopicDistribution::getAggregatedIsotopicVariants(int charge /* = 0 */)
	{
		const std::map<std::string, int> ele_count = _compo.get();
		
		std::map<std::string, int>::const_iterator it = ele_count.begin();

		// Update _var_num.
		size_t max_var = _compo.getMaxNumVariants();
		// Prevent undesirable number of variants to crash the memory.
		if(_var_num > max_var)
			_var_num = max_var;

		_iso_const.setOrder(_var_num);

		// Build the constants for all elements.
		for(; it != ele_count.end(); it++)
			_iso_const.addElement(it->first);
		
		EleSymPolyVec prob_vec = this->calculateProbabilityVec();
		std::vector<double> center_mass_vec = this->calculateCenterMassVec(prob_vec);

		AggregatedIsotopicVariants peakset;

		// Get the maximum value from the prob list.
		// TBD: For large molecules, this will be a problem, since the prob is too small.
		//double max_prob = 1e-12;
		//for(size_t i=0; i<= _var_num; i++)
		//{
		//	if(prob_vec.at(i) > max_prob)
		//		max_prob = prob_vec.at(i);
		//}

		for(size_t i=0; i<= _var_num; i++)
		{
			//if(prob_vec.at(i)/max_prob < 1e-3)
			//	break;

			double adjusted_mz = (charge == 0 ? center_mass_vec.at(i): calculateMZ(center_mass_vec.at(i), charge));

			PeakPtr peak = boost::make_shared<Peak>(adjusted_mz, prob_vec.at(i));
			peakset.addPeak(peak);
			
		}

		return peakset;
		
	}

	void IsotopicDistribution::setVaraibles( Composition& compo, int order /*= -1*/ )
	{
		setComposition(compo);
		setOrder(order);
	}

	double calculateMass( double mz, int charge)
	{	
		// Be careful of the ion mode.
		param::Param& param = param::Param::Instance();
		double electron_mass = param.getParameter<double>("electron_mass").first;

		return abs(charge) * mz - charge * (Composition("H").getMass() -electron_mass);
	}

	double calculateMZ(double mass, int charge, int pre_charge /* = 0 */)
	{
		param::Param& param = param::Param::Instance();
		double electron_mass = param.getParameter<double>("electron_mass").first;

		int coef_h = (pre_charge == 0 ? charge : pre_charge);

		// TBD: this should be controlled by parameter.
		// int coef_e = (pre_charge == 0 ? 1 : 0);
		return (mass + coef_h * (Composition("H").getMass() - electron_mass))/abs(charge);
	}
}