#include <iostream>
#include "nr3.h"
#include "ludcmp.h"
#include "utilities.h"
#include <iostream>

using namespace std;

void LUdecomposition(MatDoub &A)
{
	MatDoub L(3, 3);
	MatDoub U(3, 3);

	L[0][0] = 1;
	L[1][1] = 1;
	L[2][2] = 1;

	for (int j = 0; j < 3; j++)
	{
		for (int i = 0; i <= j; i++)
		{
			U[i][j] = A[i][j];
			for (int k = 0; k < i - 1; k++)
			{
				U[i][j] -= L[i][k] * U[k][j];
			}
		}
		for (int i = j + 1; i < 3; i++)
		{
			L[i][j] = A[i][j];
			for (int k = 0; k < j - 1; k++)
			{
				L[i][j] -= L[i][k] * U[k][j];
			}
			L[i][j] /= U[j][j];
		}
	}

	util::print(U);
	util::print(L);
	
}

int main()
{

	// Exercise 1:
	// Solve A x = b using LU decomposition, and print the result.

	MatDoub A(3, 3);
	A[0][0] = 1.0;
	A[0][1] = 0;
	A[0][2] = 0;

	A[1][0] = 0;
	A[1][1] = 1;
	A[1][2] = 0;

	A[2][0] = 0;
	A[2][1] = 0;
	A[2][2] = 1;

	LUdecomposition(A);

	VecDoub b(3);
	b[0] = 5.0;
	b[1] = 18.0;
	b[2] = 6.0;

	VecDoub x(3);
	util::print(b);
	LUdcmp ldcmp=LUdcmp(A);
	// evaluate x
	ldcmp.solve(b,x);
	// print x
	util::print(x);
	// Verify that A == LU

	return 0;
}
