#include <iostream>
#include <fstream>
#include "nr3.h"
#include "banded.h"

using namespace std;

std::ostream& operator<<(std::ostream& os, const MatDoub& x)
{
    os << "Matrix:" << "\n";
    for (size_t i = 0; i < x.nrows(); i++)
    {
        for (size_t j = 0; j < x.ncols(); j++)
        {
            os << x[i][j] << ", ";
        }
        os << "\n";
    }

    return os;
}


// General poisson equation:
// -(d^2u/d^2x + d^2u/d^2y) + lambda u(x, y) = f(x, y)

// Problem at hand:
// -(d^2u/d^2x + d^2u/d^2y) = 1 + x + y
const Doub lambda = 0;

Doub f(Doub x, Doub y) {
    return 1 + x + y;
}

template<class T>
void print(const T &v) {
    cout << setw(14) << v;
}

int main() {
    Doub x0 = 0, x1 = 1;
    Doub y0 = x0, y1 = x1;

    cout << "\numi = u(0.5, 0.5)_i\n"
            "using k=2 for richardson extrapolation error estimate\n";

    for (auto &v : {"i", "N", "umi", "umi - umi-1", "k", "e", "\n"})
        print(v);

    // write the solutions to a file for plotting
    ofstream fs("../data.csv");

    int N = 4;
    Doub umisub1 = NAN, umisub2 = NAN;
    for (int i = 0; N <= 64; i++) {
        Doub h = (x1 - x0) / N;  // same h for x and y
        // since u(x, y) is specified at the boundary (u=0), we have (N-1)^2 unknowns
        int n = N - 1;
        // A * w = phi
        // because the unknowns are coupled to their neighbours in a (N-1) x (N-1) grid,
        // the matrix becomes band-diagonal with m1 = m2 = N - 1 (see NR p. 58-61)
        // We thus build a compact band-diagonal matrix
        MatDoub A(n * n, n + 1 + n, 0.);  // rows: number of unknowns as usual, cols: total band width (m1 + 1 + m2)
        VecDoub phi(n * n);

        // for every (col, row) in the grid of unknowns
        for (int col = 0; col < n; col++) {
            Doub x = x0 + (col + 1) * h;
            for (int row = 0; row < n; row++) {
                Doub y = y0 + (row + 1) * h;

                int j = row * n + col;  // the index of the unknown
                // note that A[a][b] in the band-diagonal matrix is A[a][a + b - n] in the full matrix
                A[j][n] = 4 + h * h * lambda;  // A[j][n] thus corresponds to Ajj
                phi[j] = h * h * f(x, y);
                if (col > 0)  // interior point to the left
                    A[j][n - 1] = -1;
                if (col < n - 1)  // interior point to the right
                    A[j][n + 1] = -1;
                if (row > 0)  // interior point below
                    A[j][0] = -1;
                if (row < n - 1)  // interior point above
                    A[j][2 * n] = -1;
                // if u(x, y) != 0 at the boundary, we would have to add the appropriate boundary values to phi
            }
        }

        if (N == 4) {
            std::cout << A << std::endl;
        }

        // O(n^4) with Bandec, compared to O(n^6) with lu
        Bandec bandec(A, n, n);  // the compact band-diagonal matrix, m1, m2 (m1 = m2 = n)
        VecDoub w(n * n);
        bandec.solve(phi, w);
        Doub umi = w[(n * n) / 2];
        Doub k_est = log2((umisub2 - umisub1) / (umisub1 - umi));
        Doub e = (umisub1 - umi) / (4 - 1);

        for (auto &v : {i + 1, N})
            print(v);
        for (auto &v : {umi, umi - umisub1, k_est, e})
            print(v);
        cout << endl;

        umisub2 = umisub1;
        umisub1 = umi;
        N *= 2;

        // write the solutions to a file for plotting
        for (int j = 0; j < w.size(); j++)
            fs << w[j] << ",";
        fs << endl;
    }

    return 0;
}