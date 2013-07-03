/********************************************************************
	created:	2012/11/13
	created:	13:11:2012   9:12
	filename: 	PeakListTest.cpp
	file path:	GAG\test
	file base:	PeakListTest
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include "PeakList.h"
#include <boost/make_shared.hpp>

int main(int argc, char* argv[])
{
	using namespace std;
	using namespace brain;

	RawList peaklist;

	PeakPtr pk9 = boost::make_shared<Peak>(471.2, 1320.5);
	peaklist.addPeak(pk9);

	PeakPtr pk10 = boost::make_shared<Peak>(472.2, 170.6);
	peaklist.addPeak(pk10);

	PeakPtr pk11 = boost::make_shared<Peak>(473.2, 222222.7);
	peaklist.addPeak(pk11);

	PeakPtr pk12 = boost::make_shared<Peak>(474.2, 26.4);
	peaklist.addPeak(pk12);

	PeakPtr pk13 = boost::make_shared<Peak>(475.2, 32.5);
	peaklist.addPeak(pk13);

	// 2. printPeakList with given order.

	cout << "Sort by MZ" << endl;
	peaklist.printPeakList<peak_mz>();

	cout << "Sort by Intensity" << endl;
	peaklist.printPeakList<peak_intensity>();

	cout << "Base Peak: " << endl;
	cout << peaklist.getBasePeak<peak_intensity>() << endl;

	cout << "A+1: " << endl;
	cout << peaklist.getPeakByShift<peak_intensity>(1) << endl;

	cout << "A+2: " << endl;
	cout << peaklist.getPeakByShift<peak_intensity>(2) << endl;

	cout << "A+3: " << endl;
	cout << peaklist.getPeakByShift<peak_intensity>(3) << endl;

	cout << "A-1: " << endl;
	cout << peaklist.getPeakByShift<peak_intensity>(-1) << endl;

	cout << "A-2: " << endl;
	cout << peaklist.getPeakByShift<peak_intensity>(-2) << endl;

	cout << "A-3: " << endl;
	cout << peaklist.getPeakByShift<peak_intensity>(-3) << endl;

	cout << "The End!" << endl;
	cin.get();

}