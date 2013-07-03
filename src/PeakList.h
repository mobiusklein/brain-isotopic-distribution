/********************************************************************
	created:	2012/08/03
	created:	3:8:2012   15:21
	filename: 	PeakList.h
	file path:	GAGPL\SPECTRUM
	file base:	PeakList
	file ext:	h
	author:		Han Hu
	
	purpose:	
*********************************************************************/
#ifndef BRAIN_PEAKLIST_H
#define BRAIN_PEAKLIST_H

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/typeof/std/utility.hpp>
#include "Peak.h"

#include <iostream>
#include <cmath>
//#include "dataanalysis.h"

namespace brain
{
	using namespace ::boost;
	using namespace ::boost::multi_index;

	// Create tags.
	struct peak_mz{};
	struct peak_intensity{};

	typedef 
	multi_index_container<
		PeakPtr,
		indexed_by<
			ordered_unique<
				tag<peak_mz>, member<Peak, double, &Peak::mz> >,
			ordered_non_unique<
				tag<peak_intensity>, member<Peak, double, &Peak::intensity>, std::greater<double> >
		>
	> PeakContainer;


	struct change_type
	{
		change_type(const std::string& new_type)
			:new_type(new_type)	
		{}

		void operator()(std::string& str)
		{
			str = new_type;
		}

	private:
		std::string new_type;
	};

	template<typename PeakType, typename MultiIndexContainer>
	class PeakList
	{

	public:
		inline void addPeak(const shared_ptr<PeakType> pk)
		{
			_peakset.insert(pk);
		}
		inline size_t getSize() const
		{
			return _peakset.size();
		}

		template<typename AbundanceTag>
		inline boost::shared_ptr<PeakType> getBasePeak()
		{
			return getPeakByShift<AbundanceTag>(0);
		}

		inline boost::shared_ptr<PeakType> getMonoisotopicPeak()
		{
			return *(_peakset.get<peak_mz>().begin());
		}

		// Find the peak with maximum value, and follow the mz axis to get upstream/downstream peaks. The shift is important for theoretical peak cluster.
		template<typename AbundanceTag>
		boost::shared_ptr<PeakType> getPeakByShift(const int& shift)
		{
			typedef MultiIndexContainer::index<AbundanceTag>::type peakset_by_type;
			typedef MultiIndexContainer::index<peak_mz>::type peakset_by_mz;

			// Get the base peak.
			peakset_by_type& type_index = _peakset.get<AbundanceTag>();
			peakset_by_mz& mz_index = _peakset.get<peak_mz>();
			peakset_by_type::iterator it1 = type_index.begin();
			
			// Project the iterator to mz iterator.
			if(shift == 0) {
				return *it1;
			} else {
				// Estimate if the shift is valid.
				peakset_by_mz::iterator it2 = _peakset.project<peak_mz>(it1);
				int t = shift > 0 ? std::distance(it2, mz_index.end())-1 : std::distance(mz_index.begin(), it2);

				if(t < abs(shift))
					return boost::make_shared<PeakType>();

				std::advance(it2, shift);
				return *it2;
			}
		}

		template<typename AbundanceTag>
		double getMassDifferenceByShift(const int& shift1, const int& shift2)
		{
			return this->getPeakByShift<AbundanceTag>(shift2)->mz - this->getPeakByShift<AbundanceTag>(shift1)->mz;
		}

		template<typename Tag>
		void print(const typename MultiIndexContainer::index<Tag>::type& i)
		{
			std::copy(i.begin(),i.end(),std::ostream_iterator<boost::shared_ptr<PeakType> >(std::cout));
		}

		// Print the peaklist based on specified order.
		template<typename Tag>
		void printPeakList(std::ostream& os = std::cout)
		{
			/* obtain a reference to the index tagged by Tag */  
			const MultiIndexContainer::index<Tag>::type& i = _peakset.get<Tag>();  

			/* dump the elements of the index to cout */
			std::copy(i.begin(),i.end(),std::ostream_iterator<boost::shared_ptr<PeakType> >(os)); 
		}

		template<typename Tag>
		typename MultiIndexContainer::index<Tag>::type& getPeakListByType()
		{
			return _peakset.get<Tag>();		
		}

		inline MultiIndexContainer& getPeakContainer()
		{
			return _peakset;
		}

	private:
		MultiIndexContainer _peakset;
	};

	typedef PeakList<Peak, PeakContainer> RawList;
	
	typedef RawList AggregatedIsotopicVariants;
	
}
#endif /* BRAIN_PEAKLIST_H */