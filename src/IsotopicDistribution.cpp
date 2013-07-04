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

		// Iterate over all elements.
		for(; it != ele_count.end(); it++)
		{
			phi_value += _iso_const.getNthElementPowerSum(it->first, order) * it->second;
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
		
		NewtonGirardFormulae newton(_compo.getMaxNumVariants());
		EleSymPolyVec prob_vec;
		newton.updateParameters(phi_vec, prob_vec);
		//newton.setPowerSum(phi_vec);

		//EleSymPolyVec prob_vec = newton.getElementarySymmetricPolyVec();
		for(size_t i=0; i<prob_vec.size(); i++)
		{
			int sign = (i%2 == 0 ? 1 : -1);
			// q(j) = q(0) * e(j) * (-1)^j.
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
			// For the specified symbol, the value of the coefficient will reduce 1.
			int coef = (it->first == symbol ? it->second - 1 : it->second);
			phi_value += _iso_const.getNthElementPowerSum(it->first, order) * coef;
		}

		// The part of v(r^-l) now becomes (v-1)(r^-l) + r'^-1.
		phi_value += _iso_const.getNthModifiedElementPowerSum(symbol, order); 

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

		std::map<std::string, int> ele_count = _compo.get();

		PeriodicTable& ptable = PeriodicTable::Instance();

		NewtonGirardFormulae newton(_compo.getMaxNumVariants());

		std::map<std::string, EleSymPolyVec> esp_map;
		for(std::map<std::string, int>::const_iterator it = ele_count.begin();
			it != ele_count.end(); it++)
		{
			const Isotope& mono = ptable.getIsotopeByRelativeShift(it->first);
			PowerSumVec& ps_vec = this->calculateModifiedPhiValueVec(it->first);
			EleSymPolyVec esp_vec;
			newton.updateParameters(ps_vec, esp_vec);
			
			esp_map.insert(std::make_pair(it->first, esp_vec));
		}

		for(size_t i = 0; i <= _var_num; i++)
		{
			std::map<std::string, EleSymPolyVec>::iterator it = esp_map.begin();
			int sign = (i%2 == 0 ? 1 : -1);
			double center_mass = 0.0;
			for(; it != esp_map.end(); it++)
			{
				const Isotope& mono = ptable.getIsotopeByRelativeShift(it->first);

				center_mass += ele_count[it->first] * sign * it->second.at(i) * monopeak.intensity * mono.mass;
			}
			// m(j) = sum(m(jk) * p(jk))/sum(p(jk))
			mass_vec.push_back(center_mass/prob_vec.at(i));
		}

		//for(size_t i=0; i<=_var_num; i++)
		//{
		//	std::map<std::string, int>::const_iterator it = ele_count.begin();
		//	
		//	int sign = (i%2 == 0 ? 1 : -1);

		//	double center_mass = 0.0;

		//	for(; it != ele_count.end(); it++)
		//	{				
		//		const Isotope& mono = ptable.getIsotopeByRelativeShift(it->first);
		//		
		//		PowerSumVec& ps_vec = this->calculateModifiedPhiValueVec(it->first);
		//		EleSymPolyVec esp_vec;
		//		newton.updateParameters(ps_vec, esp_vec);
		//		//newton.setPowerSum(this->calculateModifiedPhiValueVec(it->first));
		//		
		//		// For each element, multiply the number of atoms by the modified phi.
		//		//center_mass += it->second * sign * newton.getElementarySymmetricPoly(i) * monopeak.intensity * mono.mass;
		//		center_mass += it->second * sign * esp_vec.at(i) * monopeak.intensity * mono.mass;

		//	}
		//	// m(j) = sum(m(jk) * p(jk))/sum(p(jk))
		//	mass_vec.push_back(center_mass/prob_vec.at(i));
		//}

		return mass_vec;
	}

	void IsotopicDistribution::setOrder(int order)
	{
		// Controlling the number of variants, and prevent unnecessary calculation.
		size_t max_num = _compo.getMaxNumVariants();
		if(order == -1)
			_var_num = max_num;
		else
			_var_num = (size_t)order > max_num ? max_num : order;
		
		// Adding elements from the composition into _iso_const.
		updateIsotopicConstants();
	}

	AggregatedIsotopicVariants IsotopicDistribution::getAggregatedIsotopicVariants(int charge /* = 0 */)
	{
		
		EleSymPolyVec prob_vec = this->calculateProbabilityVec();
		std::vector<double> center_mass_vec = this->calculateCenterMassVec(prob_vec);

		AggregatedIsotopicVariants peakset;

		avg_mass = 0.0; 
		double sum(0.0);

		for(size_t i=0; i<= _var_num; i++)
		{
			double adjusted_mz = (charge == 0 ? center_mass_vec.at(i): calculateMZ(center_mass_vec.at(i), charge));

			PeakPtr peak = boost::make_shared<Peak>(adjusted_mz, prob_vec.at(i));
			peakset.addPeak(peak);	

			// Code for calculating avg_mass
			avg_mass += adjusted_mz * prob_vec.at(i);
			sum += prob_vec.at(i);
		}

		avg_mass /= sum;

		return peakset;
		
	}

	void IsotopicDistribution::updateIsotopicConstants()
	{
		const std::map<std::string, int> ele_count = _compo.get();

		std::map<std::string, int>::const_iterator it = ele_count.begin();

		// Build the constants for all elements.
		for(; it != ele_count.end(); it++)
			_iso_const.addElement(it->first);

		_iso_const.updateOrder(_var_num);
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