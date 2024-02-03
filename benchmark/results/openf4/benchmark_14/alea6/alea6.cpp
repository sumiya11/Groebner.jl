// alea6
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

	variableName.push_back("x0");
	variableName.push_back("x1");
	variableName.push_back("x2");
	variableName.push_back("x3");
	variableName.push_back("x4");
	variableName.push_back("x5");
	polynomialArray.emplace_back("5*x0^2*x3 + 37*x1*x3*x4 + 32*x1*x3*x5 + 21*x3*x5 + 55*x4*x5");
	polynomialArray.emplace_back("23*x1^2*x4 + 57*x1*x2*x4 + 56*x1*x4^2 + 39*x0*x1*x5 + 52*x3*x4*x5 + 10*x2^2");
	polynomialArray.emplace_back("33*x0^2*x3 + 32*x1*x3^2 + 51*x1^2*x4 + 42*x0*x3*x5 + x5^3 + 51*x0^2");
	polynomialArray.emplace_back("44*x0*x3^2 + 47*x1*x4^2 + 43*x3*x4^2 + 2*x2*x4*x5 + 42*x1*x3 + 12*x2*x3");
	polynomialArray.emplace_back("49*x0^2*x2 + 11*x0*x1*x2 + 45*x1^2*x4 + 83*x0*x3*x4 + 54*x0*x3");
	polynomialArray.emplace_back("48*x0*x2*x3 + 2*x2^2*x3 + 36*x3^3 + 59*x2^2*x5 + 17*x2 + 45*x4");

	vector<string> basis = groebnerBasisF4(1073741827, 6, variableName, polynomialArray, 1, 0);

	std::cout << "The basis contains " << basis.size() << " elements." << std::endl;


	ofstream output;
	output.open("/home/demin/Groebner.jl/benchmark/results/openf4/benchmark_14/alea6//output");
	output << "x0, x1, x2, x3, x4, x5" << endl;
	output << "1073741827" << endl;
	for (size_t i = 0; i < basis.size(); i++) {
		output << basis[i] << "," << endl;
	}
	output.close();

	return 0;
}

