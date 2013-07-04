/********************************************************************
	created:	2013/07/01
	created:	1:7:2013   11:12
	filename: 	IsotopicDistribution.h
	file path:	BRAIN\src
	file base:	IsotopicDistribution
	file ext:	h
	author:		Han Hu
	
	purpose:	BRAIN algorithm for calculating theoretical isotopic 
						distribution.
*********************************************************************/
#ifndef BRAIN_ISOTOPICDISTRIBUTION_H
#define BRAIN_ISOTOPICDISTRIBUTION_H

#include "Composition.h"
#include "IsotopicConstants.h"
#include "PeakList.h"

namespace brain
{
	// The aggregated masses and corresponding probabilities.
	// Note that the double type for probability should be converted 
	// to float type in order for the abundance matching procedure.

	// Based on the mode, charge state and mz, calculate the original mass.
	double calculateMass(double mz, int charge);

	// Calculate mz from mode, mass and charge. By default the number of proton equals to the charge state value.
	double calculateMZ(double mass, int charge, int pre_charge = 0);

	// Implementation of BRAIN algorithm.
	class IsotopicDistribution
	{
	public:
		IsotopicDistribution()
			: _iso_const(IsotopicConstants::Instance()) {}
		//IsotopicDistribution(size_t order)
		//: _iso_const(IsotopicConstants::Instance()), _var_num(order) {}
		IsotopicDistribution(Composition& compo, int order = -1)
			: _iso_const(IsotopicConstants::Instance()), _compo(compo)
		{
			// Implicitly build the distribution.
			createMonoPeak();
			setOrder(order);
		}

		EleSymPolyVec calculateProbabilityVec();

		std::vector<double> calculateCenterMassVec(const EleSymPolyVec& prob_vec);

		// The peak cluster can be adjusted by charge. Notice that the charge here can be either positive or negative.
		AggregatedIsotopicVariants getAggregatedIsotopicVariants(int signed_charge = 0);

		inline double getAverageMass() const
		{
			return avg_mass;
		}

		inline std::string getCompositionString() const
		{
			return _compo.getCompositionString();
		}

	private:
		// phi_vec[0] should be assigned 0.
		IsotopicConstants& _iso_const;

		// Composition of the molecular formula.
		Composition _compo;
		
		Peak monopeak;

		// The average mass calculated based on returned isotopic peaks. 
		// It can be used to control the number of calculated peaks.
		double avg_mass;

		// The number of isotopic variants.
		size_t _var_num;

		// The phi_vec[index] has been assigned, just return the value,
		// otherwise, calculate it from the scratch.
		double calculatePhiValue(size_t order);
		PowerSumVec calculatePhiValueVec();

		double calculateModifiedPhiValue(const std::string& symbol, size_t order);
		PowerSumVec calculateModifiedPhiValueVec(const std::string& symbol);

		// Calculating the information of monoisotopic peak. 
		void createMonoPeak();

		void setOrder(int order);

		void updateIsotopicConstants();


	};

}

#endif /* BRAIN_ISOTOPICDISTRIBUTION_H */