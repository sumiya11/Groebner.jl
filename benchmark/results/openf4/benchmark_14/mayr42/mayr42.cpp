// mayr42
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
	variableName.push_back("x33");
	variableName.push_back("x34");
	variableName.push_back("x35");
	variableName.push_back("x36");
	variableName.push_back("x37");
	variableName.push_back("x38");
	variableName.push_back("x39");
	variableName.push_back("x40");
	variableName.push_back("x41");
	variableName.push_back("x42");
	variableName.push_back("x43");
	variableName.push_back("x44");
	variableName.push_back("x45");
	variableName.push_back("x46");
	variableName.push_back("x47");
	variableName.push_back("x48");
	variableName.push_back("x49");
	variableName.push_back("x50");
	variableName.push_back("x51");
	polynomialArray.emplace_back("x4*x49 + 1073741826*x10*x51");
	polynomialArray.emplace_back("x3*x48 + 1073741826*x9*x51");
	polynomialArray.emplace_back("x2*x47 + 1073741826*x8*x51");
	polynomialArray.emplace_back("x1*x46 + 1073741826*x7*x51");
	polynomialArray.emplace_back("x4*x44 + 1073741826*x9*x49");
	polynomialArray.emplace_back("x3*x43 + 1073741826*x8*x48");
	polynomialArray.emplace_back("x2*x42 + 1073741826*x7*x47");
	polynomialArray.emplace_back("x1*x41 + 1073741826*x6*x46");
	polynomialArray.emplace_back("x4*x39 + 1073741826*x9*x49");
	polynomialArray.emplace_back("x3*x38 + 1073741826*x8*x48");
	polynomialArray.emplace_back("x2*x37 + 1073741826*x7*x47");
	polynomialArray.emplace_back("x1*x36 + 1073741826*x6*x46");
	polynomialArray.emplace_back("x9*x34 + 1073741826*x9*x49");
	polynomialArray.emplace_back("x4*x34 + 1073741826*x5*x51");
	polynomialArray.emplace_back("x8*x33 + 1073741826*x8*x48");
	polynomialArray.emplace_back("x3*x33 + 1073741826*x4*x51");
	polynomialArray.emplace_back("x7*x32 + 1073741826*x7*x47");
	polynomialArray.emplace_back("x2*x32 + 1073741826*x3*x51");
	polynomialArray.emplace_back("x6*x31 + 1073741826*x6*x46");
	polynomialArray.emplace_back("x1*x31 + 1073741826*x2*x51");
	polynomialArray.emplace_back("x9*x14*x39 + 1073741826*x9*x29*x44");
	polynomialArray.emplace_back("x8*x13*x38 + 1073741826*x8*x28*x43");
	polynomialArray.emplace_back("x7*x12*x37 + 1073741826*x7*x27*x42");
	polynomialArray.emplace_back("x6*x11*x36 + 1073741826*x6*x26*x41");
	polynomialArray.emplace_back("x6*x26^2*x46 + 1073741826*x7*x51^3");
	polynomialArray.emplace_back("x6*x11^2*x46 + 1073741826*x2*x51^3");
	polynomialArray.emplace_back("x6*x21^2*x41 + 1073741826*x6*x46*x51^2");
	polynomialArray.emplace_back("x6*x16^2*x36 + 1073741826*x6*x46*x51^2");
	polynomialArray.emplace_back("x9*x24*x30*x39*x50 + 1073741826*x9*x29*x44*x50*x51");
	polynomialArray.emplace_back("x8*x23*x29*x38*x49 + 1073741826*x8*x28*x43*x49*x51");
	polynomialArray.emplace_back("x7*x22*x28*x37*x48 + 1073741826*x7*x27*x42*x48*x51");
	polynomialArray.emplace_back("x6*x21*x27*x36*x47 + 1073741826*x6*x26*x41*x47*x51");
	polynomialArray.emplace_back("x9*x24*x25*x39*x45 + 1073741826*x9*x29*x44*x45*x51");
	polynomialArray.emplace_back("x8*x23*x24*x38*x44 + 1073741826*x8*x28*x43*x44*x51");
	polynomialArray.emplace_back("x7*x22*x23*x37*x43 + 1073741826*x7*x27*x42*x43*x51");
	polynomialArray.emplace_back("x6*x21*x22*x36*x42 + 1073741826*x6*x26*x41*x42*x51");
	polynomialArray.emplace_back("x9*x20*x24*x39*x40 + 1073741826*x9*x29*x40*x44*x51");
	polynomialArray.emplace_back("x8*x19*x23*x38*x39 + 1073741826*x8*x28*x39*x43*x51");
	polynomialArray.emplace_back("x9*x15*x24*x35*x39 + 1073741826*x9*x29*x35*x44*x51");
	polynomialArray.emplace_back("x7*x18*x22*x37*x38 + 1073741826*x7*x27*x38*x42*x51");
	polynomialArray.emplace_back("x8*x14*x23*x34*x38 + 1073741826*x8*x28*x34*x43*x51");
	polynomialArray.emplace_back("x6*x17*x21*x36*x37 + 1073741826*x6*x26*x37*x41*x51");
	polynomialArray.emplace_back("x7*x13*x22*x33*x37 + 1073741826*x7*x27*x33*x42*x51");
	polynomialArray.emplace_back("x6*x12*x21*x32*x36 + 1073741826*x6*x26*x32*x41*x51");

	vector<string> basis = groebnerBasisF4(1073741827, 51, variableName, polynomialArray, 1, 0);

	std::cout << "The basis contains " << basis.size() << " elements." << std::endl;


	ofstream output;
	output.open("/home/demin/Groebner.jl/benchmark/results/openf4/benchmark_14/mayr42//output");
	output << "x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51" << endl;
	output << "1073741827" << endl;
	for (size_t i = 0; i < basis.size(); i++) {
		output << basis[i] << "," << endl;
	}
	output.close();

	return 0;
}

