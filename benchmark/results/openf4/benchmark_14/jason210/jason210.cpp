// jason210
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <libopenf4.h>
using namespace std;

int main (int argc, char **argv)
{

	vector<string> polynomialArray;
	vector<string> variableName;

	variableName.push_back("x1");
	variableName.push_back("x2");
	variableName.push_back("x3");
	variableName.push_back("x4");
	variableName.push_back("x5");
	variableName.push_back("x6");
	variableName.push_back("x7");
	variableName.push_back("x8");
	polynomialArray.emplace_back("x1^2*x3^4 + x2^2*x4^4 + x1*x2*x3^2*x5^2 + x1*x2*x4^2*x6^2 + x1*x2*x3*x4*x5*x7 + x1*x2*x3*x4*x6*x8");
	polynomialArray.emplace_back("x2^6");
	polynomialArray.emplace_back("x1^6");

	vector<string> basis = groebnerBasisF4(1073741827, 8, variableName, polynomialArray, 1, 0);

	std::cout << "The basis contains " << basis.size() << " elements." << std::endl;


	ofstream output;
	output.open("/home/demin/Groebner.jl/benchmark/results/openf4/benchmark_14/jason210//output");
	output << "x1, x2, x3, x4, x5, x6, x7, x8" << endl;
	output << "1073741827" << endl;
	for (size_t i = 0; i < basis.size(); i++) {
		output << basis[i] << "," << endl;
	}
	output.close();

	return 0;
}

