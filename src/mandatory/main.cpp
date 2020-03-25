#include <iostream>
#include "nr3.h"
#include "utilities.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <svd.h>


constexpr auto NUM_ROWS = 500;
constexpr auto NUM_COLS = 4;

double l2(const VecDoub & x)
{
    double sum = 0;
    for (int i = 0; i < x.size(); i++)
        sum += std::pow(x[i],2);
    return sqrt(sum);
}

VecDoub operator-(const VecDoub& a, const VecDoub& b)
{
    int len = a.size() > b.size() ? a.size() : b.size();
    VecDoub c(len);
    for (int i = 0; i < len; i++)
        c[i] = a[i] - b[i];
    return c;
}

VecDoub std_dev_q(const SVD& svd, int tsh_num = 0)
{
    VecDoub d_x(svd.w.size());
    for (int i = 0; i < d_x.size(); i++)
    {
        double sum = 0;
        for (int j = 0; j < svd.w.size() - tsh_num; j++)
            sum += std::pow(svd.v[i][j] / (svd.w[j]), 2);
        d_x[i] = sqrt(sum);
    }
    return d_x;
}

double residual_error(const MatDoub& A, VecDoub& q, VecDoub& z)
{
    return (l2 (A * q - z) / l2(z));
}

void generate_A(const MatDoub& data, MatDoub& A, const double sigma = 1)
{
    for (int i = 0; i < data.nrows(); i++)
    {
        A[i * 2][0]     = 1/sigma;
        A[i * 2][1]     = 0/sigma;
        A[i * 2][2]     = std::cos(data[i][0])/sigma;
        A[i * 2][3]     = std::cos(data[i][0] + data[i][1])/sigma;

        A[i * 2 + 1][0] = 0/sigma;
        A[i * 2 + 1][1] = 1/sigma;
        A[i * 2 + 1][2] = std::sin(data[i][0])/sigma;
        A[i * 2 + 1][3] = std::sin(data[i][0] + data[i][1])/sigma;
    }
}

void generate_b(const MatDoub &data, VecDoub& b, const double sigma = 1)
{
    for (int i = 0; i < data.nrows(); i++)
    {
        b[i * 2]     = data[i][2]/sigma;
        b[i * 2 + 1] = data[i][3]/sigma;
    }
}

int main(int argc, char const *argv[])
{
    MatDoub data1(NUM_ROWS,NUM_COLS);
    MatDoub data2(NUM_ROWS,NUM_COLS);

    std::cout << "Please make sure to give the correct path to d1 and d2 <3" << std::endl;
    ifstream d1("/home/daniel/Desktop/numerical-methods/data/d1");
    ifstream d2("/home/daniel/Desktop/numerical-methods/data/d2");

    for (int i = 0; i < NUM_ROWS; i++)
    {
        for (int j = 0; j < NUM_COLS; j++)
        {
            d1 >> data1[i][j];
            d2 >> data2[i][j];
        }
    }

    d1.close(); d2.close();

    MatDoub A1(2 * NUM_ROWS, NUM_COLS);
    MatDoub A2(2 * NUM_ROWS, NUM_COLS);
    VecDoub z1(2 * NUM_ROWS);
    VecDoub z2(2 * NUM_ROWS);
    
    double sigma = 0.1;

    generate_A(data1,A1,sigma);
    generate_b(data1,z1,sigma);

    generate_A(data2,A2,sigma);
    generate_b(data2,z2,sigma);

    SVD svd1(A1);
    SVD svd2(A2);

    std::cout << "Machine precision := " << svd1.eps << std::endl;

    util::print(svd1.w,"w1");
    std::cout << "The condition of the design matrix A1: " << svd1.inv_condition() << std::endl;
    std::cout << "\n";

    util::print(svd2.w,"w2");
    std::cout << "The condition of the design matrix A2: " << svd2.inv_condition() << std::endl;
    std::cout << "\n";

    VecDoub q1(A1.ncols());
    VecDoub q2(A2.ncols());
    VecDoub e_residual2(svd2.w.size());

    svd1.solve(z1, q1);
    util::print(q1, "q1");
    std::cout << "Residual error 1: " << residual_error(A1,q1,z1) << std::endl;
    std::cout << "\n";
    util::print(std_dev_q(svd1), "std_dev_q(scd1)");
    std::cout << "\n";

    svd2.solve(z2,q2);
    util::print(q2,"q2");
    std::cout << "Residual error 2 (without threshold): " << residual_error(A2,q2,z2) << std::endl;
    std::cout << "\n";
    util::print(std_dev_q(svd2), "std_dev_q(svd2)");
    std::cout << "\n";

    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";

    double THRESHOLD_LIMIT = svd2.w[0];
    int i = 0;

    for (double threshhold = svd2.w[svd2.w.size() - 1]; threshhold < THRESHOLD_LIMIT; threshhold = svd2.w[svd2.w.size() - 1 - i])
    {
        svd2.solve(z2, q2, svd2.w[3]);
        util::print(q2, "q2");
        e_residual2[i] = residual_error(A2,q2,z2);
        i++;
        std::cout << "\n";
        util::print(std_dev_q(svd2,i), "std_dev_q(svd2)");
        std::cout << "\n";
    }

    util::print(e_residual2, "Residual error 2: ");

    return 0;
}
