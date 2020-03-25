#include <iostream>
#include "nr3.h"
#include "ludcmp.h"
#include "cholesky.h"
#include "utilities.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <svd.h>
#include <exception>
#include <eigen_sym.h>

#define GRAM 1

void calcParam(const std::vector<double_t> & Y, const std::vector<double_t> & X, MatDoub & A, MatDoub & b, int o);
void svdProblem(const MatDoub & A, const MatDoub &b);

int main() 
{

	std::vector<double_t> Y;
	std::vector<double_t> X;
    std::ifstream inFile;
	double_t data;
    MatDoub A;
    MatDoub b;

    // Open file
    inFile.open("/home/daniel/Desktop/numerical-methods/data/FilipData.dat");

    // Is file able to be opened?
    if (!inFile) 
	{
        cout << "Unable to open file";
        exit(1); // terminate with error
    }

    while (inFile >> data) 
	{
		Y.push_back(data);
		inFile >> data;
		X.push_back(data);
    }

    inFile.close();
    calcParam(Y, X, A, b, 10);
    svdProblem(A,b);
	
    return 0;
}

void calcParam(const std::vector<double_t> & Y, const std::vector<double_t> & X, MatDoub & A, MatDoub & b, int o)
{
    A.resize(Y.size(),o);
    for (size_t i = 0; i < A.nrows(); i++)
        for (size_t j = 0; j < A.ncols(); j++)
            A[i][j] = std::pow(X[i],j);

    b.resize(Y.size(),1);
    for (size_t i = 0; i < b.nrows(); i++)
        b[i][0] = Y[i]; 
}



void svdProblem(const MatDoub & A, const MatDoub &b)
{


    printf("A: The size is %d %d \n", A.nrows(), A.ncols());
    printf("b: The size is %d %d \n", b.nrows(), b.ncols());

    MatDoub x(A.ncols(),1);

    // Solution by SVD
    SVD svdcmp(A);
    std::cout << "Solution by SVD" << std::endl;
    svdcmp.solve(b, x, -1);
    util::print(x);
    util::print(svdcmp.w);
}
