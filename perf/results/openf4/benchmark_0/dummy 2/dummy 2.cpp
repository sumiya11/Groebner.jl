// dummy 2
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
	polynomialArray.emplace_back("x^2 + y^2 + 6");
	polynomialArray.emplace_back("x + y");

	vector<string> basis = groebnerBasisF4(7, 2, variableName, polynomialArray, 1, 0);

	std::cout << "The basis contains " << basis.size() << " elements." << std::endl;

	return 0;
}

