#include <iostream>
#include <fstream>
#include <functional>
#include "nr3.h"
#include "tridag.h"

using namespace std;

struct Config
{
    string name;
    Doub atol = 1e-4;
    Doub x0 = 0, x1 = 1, t0 = 0, t1 = 20;
    Doub alpha = 1; // the diffusion constant, not the richardson alpha.
    function<Doub(Doub)> a = [](Doub t) { return 0; };
    function<Doub(Doub)> b = [](Doub t) { return 1; };
    function<Doub(Doub)> g = [](Doub x) { return pow(x, 4); };
    function<Doub(Doub, Doub)> f = [](Doub x, Doub t) { return 10 * x * (1 - x) * exp(-t / 10); };
};

template <class T>
void print(const T &v)
{
    cout << setw(14) << v;
}

// for plotting
void write_vec(ofstream &of, const VecDoub &vec)
{
    for (int i = 0; i < vec.size(); i++)
        of << vec[i] << ",";
    of << endl;
}

std::ostream &operator<<(std::ostream &os, const MatDoub &x)
{
    os << "Matrix:"
       << "\n";
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

std::ostream &operator<<(std::ostream &os, const VecDoub &x)
{
    os << "Vector:"
       << "\n";

    for (size_t i = 0; i < x.size(); i++)
    {
        os << x[i] << ", ";
    }

    return os << "\n";
}

int main()
{
    Config c_default = Config();
    c_default.name = "default (alpha=1, b=1, f=10x(1-x)exp(-t/10))";
    // the large diffusion constant, alpha, compared to x1-x0 and t1-t0,
    // makes it a bit hard to see the diffusion in the plots.
    // the below problems are chosen to display more intuitive diffusion examples:
    Config c_smooth = Config();
    c_smooth.name = "smooth (alpha=1e-2, b=1, f=0)";
    c_smooth.alpha = 1e-2;
    c_smooth.f = [](Doub x, Doub t) { return 0; };

    Config c_periodic = Config();
    c_periodic.name = "periodic (alpha=1e-2, b=.5+.5*cos(t), f=0)";
    c_periodic.alpha = c_smooth.alpha;
    c_periodic.f = c_smooth.f;
    c_periodic.b = [](Doub t) { return .5 + .5 * cos(t); };

    vector<Config> configs = {c_default, c_smooth, c_periodic};

    ofstream fs("../data.csv"); // for plotting

    cout << "\numi = u(0.5, 20)_i\n"
            "using k=2 for richardson extrapolation error estimate\n\n";

    for (auto &c : configs)
    {

        fs << c.name << endl;
        fs << c.x0 << "," << c.x1 << "," << c.t0 << "," << c.t1 << endl;
        cout << c.name << endl;
        for (auto &v : {"i", "Nx", "umi", "umi - umi-1", "k", "e", "\n"})
            print(v);

        int Nx = 2;
        Doub e = NAN;
        Doub umisub1 = NAN, umisub2 = NAN;
        for (int i = 0; isnan(e) or abs(e) > c.atol; i++)
        {
            Doub dx, dt;
            dt = dx = (c.x1 - c.x0) / Nx; // c = 1
            Doub r = c.alpha * dt / (dx * dx);
            int Nt = (int)round((c.t1 - c.t0) / dt); // number of steps (in time)
            // note that we can only 'hit' t1's which are multiples of dt
            // setting dt = dx is okay here, since (t1 - t0) % (x1 - x0) = 0 for our x0, x1, t0, t1

            // Nx - 1 unknowns per step + 2 for the boundaries
            // (we can let the boundaries be part of the matrix to simplify the loop later.
            //  we could also solve it keeping the boundaries out of u,
            //  but then the loop becomes a bit more complicated)
            int n = Nx + 1;
            VecDoub u(n);
            // initialize u with g(x)
            for (int j = 0; j < n; j++)
                u[j] = c.g(c.x0 + j * dx);

            fs << Nt + 1 << endl; // for plotting

            // do Nt steps
            for (int k = 0; k < Nt; k++)
            {
                write_vec(fs, u); // for plotting
                Doub tk = c.t0 + k * dt;
                Doub tkp1 = tk + dt;
                // because the unknowns are only coupled to their neighbours in one dimension,
                // the matrix is tri-diagonal and we can solve it in O(N) time
                VecDoub at(n), bt(n), ct(n), rt(n);

                // boundaries
                bt[0] = bt[n - 1] = 1;
                ct[0] = at[n - 1] = 0;
                rt[0] = c.a(tkp1);
                rt[n - 1] = c.b(tkp1);

                // the unknowns
                for (int j = 1; j < n - 1; j++)
                {
                    Doub xj = c.x0 + j * dx;
                    bt[j] = 1 + r;
                    at[j] = ct[j] = -.5 * r;
                    rt[j] = (1 - r) * u[j] +
                            .5 * r * (u[j - 1] + u[j + 1]) +
                            .5 * dt * (c.f(xj, tkp1) + c.f(xj, tk));
                }
                tridag(at, bt, ct, rt, u);
            }
            write_vec(fs, u); // for plotting

            Doub umi = u[n / 2];
            Doub k_est = log2((umisub2 - umisub1) / (umisub1 - umi));
            e = (umisub1 - umi) / (4 - 1);

            for (auto &v : {i + 1, Nx})
                print(v);
            for (auto &v : {umi, umi - umisub1, k_est, e})
                print(v);
            cout << endl;

            umisub2 = umisub1;
            umisub1 = umi;
            Nx *= 2;
        }
        fs << "end" << endl; // for plotting
    }

    cout << "\nrun plot.py for visualization\n";

    return 0;
}