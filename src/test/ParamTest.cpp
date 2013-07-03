#include "Param.h"
#include <iostream>

using namespace param;

void printParam(Param& param)
{
	// Print parameters.
	int max_charge = param.getParameter<int>("max_charge").first;
	int min_charge = param.getParameter<int>("min_charge").first;
	double tol = param.getParameter<double>("tolerance").first;
	int notexist = param.getParameter<int>("nothing").first;

	std::cout << "Max: " << max_charge << std::endl;
	std::cout << "Min: " << min_charge << std::endl;
	std::cout << "Tolerance: " << tol << std::endl;
	std::cout << "Nothing: " << notexist << std::endl;	
}

int main(int argc, char* argv[])
{
	// Test 1: Load default parameters.
	Param& param = Param::Instance();
	std::cout << "Get default value: " << std::endl;
	printParam(param);

	// Test 2: set parameters.
	std::cout << "Set parameters: " << std::endl;
	param.setParameter("max_charge", 3);
	param.setParameter("min_charge", 1);
	param.setParameter<double>("tolerance", 5e-6);
	param.setParameter<double>("A2tolerance", 1e-5);

	printParam(param);

	// Test 3: Change parameters.
	std::cout << "Change parameters:" << std::endl;
	param.setParameter("max_charge", 4);
	printParam(param);

	std::cin.get();
}