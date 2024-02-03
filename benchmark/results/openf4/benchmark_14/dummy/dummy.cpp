// dummy
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

	variableName.push_back("x");
	variableName.push_back("y");
	polynomialArray.emplace_back("x^2 + y^2 + 1073741826");
	polynomialArray.emplace_back("x + y");

	vector<string> basis = groebnerBasisF4(1073741827, 2, variableName, polynomialArray, 1, 0);

	std::cout << "The basis contains " << basis.size() << " elements." << std::endl;


	ofstream output;
	output.open("/home/demin/Groebner.jl/benchmark/results/openf4/benchmark_14/dummy//output");
	output << "x, y" << endl;
	output << "1073741827" << endl;
	for (size_t i = 0; i < basis.size(); i++) {
		output << basis[i] << "," << endl;
	}
	output.close();

	return 0;
}

