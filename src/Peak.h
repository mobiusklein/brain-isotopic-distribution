/********************************************************************
	created:	2013/07/01
	created:	1:7:2013   14:01
	filename: 	Peak.h
	file path:	BRAIN\src
	file base:	Peak
	file ext:	h
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#ifndef BRAIN_PEAK_H
#define BRAIN_PEAK_H

#include <boost/shared_ptr.hpp>
#include <iostream>

namespace brain
{
	// Forward declaration.
	struct Peak;
	typedef boost::shared_ptr<Peak> PeakPtr;
	
	struct Peak
	{
		double mz;
		double intensity;

		Peak(const double& pk_mz, const double& pk_intensity)
			: mz(pk_mz), intensity(pk_intensity)
		{}

		Peak()
			: mz(0.0), intensity(0.0)
		{}

		// Sort peak by m/z value or intensity.
		static inline bool mzSmaller(const Peak& left, const Peak& right)
		{
			return (left.mz <= right.mz);
		}
		static inline bool mzLarger(const Peak& left, const Peak& right)
		{
			return (left.mz > right.mz);
		}
		static inline bool intensitySmaller(const Peak& left, const Peak& right)
		{
			return (left.intensity <= right.intensity);
		}
		static inline bool intensityLarger(const Peak& left, const Peak& right)
		{
			return (left.intensity > right.intensity);
		}
	
		virtual inline bool empty()
		{
			return mz == 0.0 && intensity == 0.0;
		}

		virtual void printout(std::ostream& os) const;

		friend std::ostream& operator<<(std::ostream& os, const PeakPtr& pk)
		{
			pk->printout(os);
			return os;
		}
	};

}

#endif /* BRAIN_PEAK_H */