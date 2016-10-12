// Attempting to merged mutation and recombination into
// a single step

#include <algorithm>
#include <chrono>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <limits>
#include <vector>
#include "mix.h"

using namespace std;

int
main(int argc, char **argv)
{
    int argn = 1;
    unsigned seed = static_cast<unsigned>(atoi(argv[argn++]));
    double mg1 = static_cast<double>(atof(argv[argn++]));
    double mg2 = static_cast<double>(atof(argv[argn++]));
    double mnm = static_cast<double>(atof(argv[argn++]));
    double mbrk = static_cast<double>(atof(argv[argn++]));

    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed);
    vector<vector<unsigned>> g1s, g2s, mutss, brks;
    for (unsigned i = 0; i < 5000; ++i)
        {
            auto g1 = unique_fill(r, gsl_ran_poisson(r, mg1));
            auto g2 = unique_fill(r, gsl_ran_poisson(r, mg2));
            auto muts = unique_fill(r, gsl_ran_poisson(r, mnm));
            auto brk = unique_fill(r, gsl_ran_poisson(r, mbrk));
            brk.push_back(numeric_limits<unsigned>::max());
            auto output1 = mutrec_current(g1, g2, muts, brk);
            auto output2 = mutrec_new(g1, g2, muts, brk);
            if (output1 != output2)
                {
                    cerr << "failure after " << i << " successes\n";
                }
            g1s.emplace_back(std::move(g1));
            g2s.emplace_back(std::move(g2));
            mutss.emplace_back(std::move(muts));
            brks.emplace_back(std::move(brk));
            // cout << "check: " << (output1 == output2) << '\n';
        }
}
