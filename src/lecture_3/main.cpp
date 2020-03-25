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

void gram_smith();
void calcParam(const std::vector<double_t> & Y, const std::vector<double_t> & X, MatDoub & A, MatDoub & b, int o);
void transpose(const MatDoub & A, MatDoub &A_T);
void leastSquareProblem(const MatDoub & A, const MatDoub &b);

int main() 
{
    if(true)
    {
        gram_smith();
        return 0;
    }

	std::vector<double_t> Y;
	std::vector<double_t> X;
    std::ifstream inFile;
	double_t data;
    MatDoub A;
    MatDoub b;

    // Open file
    inFile.open("/home/daniel/Desktop/numerical-methods/data/PontiusData.dat");

    // Is file able to be opened?
    if (!inFile) 
	{
        cout << "Unable to open file";
        exit(1); // terminate with error
    }

    // Data is loaded
    while (inFile >> data) 
	{
		Y.push_back(data);
		inFile >> data;
		X.push_back(data);
    }

    inFile.close();

    calcParam(Y, X, A, b, 3);
    leastSquareProblem(A,b);

    X.clear();
    Y.clear();

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
    leastSquareProblem(A,b);
	
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

void transpose(const MatDoub & A, MatDoub &A_T)
{
    A_T.resize(A.ncols(),A.nrows());
    for (size_t i = 0; i < A.nrows(); i++)
    {
        for (size_t j = 0; j < A.ncols(); j++)
        {
            A_T[j][i] = A[i][j];
        }
    }
}

void gram_smith()
{
    // some array of vectors
    std::vector<VecDoub> a;

    // some test data
    double data1[] = {1.0, 1.0, 0};
    double data2[] = {3.0, 2.0, 1};

    // push data into array
    a.push_back(VecDoub(3,data1));
    a.push_back(VecDoub(3,data2));

    //
    //util::print(a[0]);
    //util::print(a[1]);

    // This lambda should be called for each vector (data...)
    auto lambda = [](std::vector<VecDoub> & b, int it)
    { 
        // This lambda takes care of the project
        auto lambda2 = [](const VecDoub & a, const VecDoub & b, const VecDoub & c)
        {
            double dot = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
            auto g = c;
            g[0] = g[0]*dot; g[1] = g[1]*dot ; g[2] = g[2]*dot;
            return g;
        };

        // This lambda subtracts two vectors
        auto lambda3 = [](const VecDoub & a, const VecDoub & b)
        {
            VecDoub c = b;
            c[0] = a[0]-b[0];  c[1] = a[1]-b[1] ;  c[2] = a[2]-b[2];
            return c;
        };

        auto e = b[it];

        for (size_t i = 0; i < it; i++)
        {
            e = lambda3(e,lambda2(b[it],b[i],b[i]));
        }
        
        double norm = sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
        
        // ERROR!?!??
        if (norm == 0)
            ;
        else
        {
            // Normalize
            e[0] = e[0] / norm;
            e[1] = e[1] / norm;
            e[2] = e[2] / norm;
            b[it] = e;
        }
    };
    
    lambda(a,0);
    lambda(a,1);

    util::print(a[0]);
    util::print(a[1]);
}



void leastSquareProblem(const MatDoub & A, const MatDoub &b)
{
    // Transpose the Matrix
    MatDoub A_T;
    transpose(A,A_T);

    MatDoub C = A_T * A;
    MatDoub c = A_T * b;

    // Solution by Cholesky
    VecDoub_IO aa(c.nrows());
    VecDoub_IO cc(c.nrows());
    for (size_t i = 0; i < c.nrows(); i++)
        cc[i] = c[i][0];

    // Solution by LU decomposition
    std::cout << "Solution by LU decomposition" << std::endl;
    LUdcmp ldcmp= LUdcmp(C);
    ldcmp.solve(cc,aa);
    util::print(aa);

    // Solution by SVD
    SVD svdcmp(C);
    std::cout << "Solution by SVD" << std::endl;
    svdcmp.solve(cc, aa);
    util::print(aa);

    // The eigen values of C
    Symmeig eigen_vals(C);
    std::cout << "Every positive semidefinite matrix only has eigenvalues ≥ 0" << std::endl;
    std::cout << "Eigenvalues of C:";
    util::print(eigen_vals.d);
}
