// bayes148
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
	variableName.push_back("x9");
	variableName.push_back("x10");
	variableName.push_back("x11");
	variableName.push_back("x12");
	variableName.push_back("x13");
	variableName.push_back("x14");
	variableName.push_back("x15");
	variableName.push_back("x16");
	variableName.push_back("x17");
	variableName.push_back("x18");
	variableName.push_back("x19");
	variableName.push_back("x20");
	variableName.push_back("x21");
	variableName.push_back("x22");
	variableName.push_back("x23");
	variableName.push_back("x24");
	variableName.push_back("x25");
	variableName.push_back("x26");
	variableName.push_back("x27");
	variableName.push_back("x28");
	variableName.push_back("x29");
	variableName.push_back("x30");
	variableName.push_back("x31");
	variableName.push_back("x32");
	polynomialArray.emplace_back("x24*x31 + 1073741826*x23*x32");
	polynomialArray.emplace_back("x24*x30 + 1073741826*x22*x32");
	polynomialArray.emplace_back("x23*x30 + 1073741826*x22*x31");
	polynomialArray.emplace_back("x24*x29 + 1073741826*x21*x32");
	polynomialArray.emplace_back("x23*x29 + 1073741826*x21*x31");
	polynomialArray.emplace_back("x22*x29 + 1073741826*x21*x30");
	polynomialArray.emplace_back("x16*x28 + 1073741826*x12*x32");
	polynomialArray.emplace_back("x20*x27 + 1073741826*x19*x28");
	polynomialArray.emplace_back("x15*x27 + 1073741826*x11*x31");
	polynomialArray.emplace_back("x20*x26 + 1073741826*x18*x28");
	polynomialArray.emplace_back("x19*x26 + 1073741826*x18*x27");
	polynomialArray.emplace_back("x14*x26 + 1073741826*x10*x30");
	polynomialArray.emplace_back("x20*x25 + 1073741826*x17*x28");
	polynomialArray.emplace_back("x19*x25 + 1073741826*x17*x27");
	polynomialArray.emplace_back("x18*x25 + 1073741826*x17*x26");
	polynomialArray.emplace_back("x13*x25 + 1073741826*x9*x29");
	polynomialArray.emplace_back("x8*x20 + 1073741826*x4*x24");
	polynomialArray.emplace_back("x18*x19 + 1073741826*x17*x20 + 1073741826*x20*x21 + x19*x22 + x18*x23 + x22*x23 + 1073741826*x17*x24 + 1073741826*x21*x24 + 1073741826*x24*x25 + x23*x26 + 2*x18*x27 + x22*x27 + x26*x27 + 1073741825*x17*x28 + 1073741826*x21*x28 + 1073741826*x25*x28 + 1073741826*x20*x29 + 1073741826*x28*x29 + x19*x30 + x27*x30 + x18*x31 + 2*x22*x31 + x26*x31 + x30*x31 + 1073741826*x17*x32 + 1073741825*x21*x32 + 1073741826*x25*x32 + 1073741826*x29*x32");
	polynomialArray.emplace_back("x7*x19 + 1073741826*x3*x23");
	polynomialArray.emplace_back("x6*x18 + 1073741826*x2*x22");
	polynomialArray.emplace_back("x5*x17 + 1073741826*x1*x21");

	vector<string> basis = groebnerBasisF4(1073741827, 32, variableName, polynomialArray, 1, 0);

	std::cout << "The basis contains " << basis.size() << " elements." << std::endl;


	ofstream output;
	output.open("/home/demin/Groebner.jl/benchmark/results/openf4/benchmark_14/bayes148//output");
	output << "x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32" << endl;
	output << "1073741827" << endl;
	for (size_t i = 0; i < basis.size(); i++) {
		output << basis[i] << "," << endl;
	}
	output.close();

	return 0;
}

