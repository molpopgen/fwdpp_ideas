// Attempting to merged mutation and recombination into
// a single step

#include <algorithm>
#include <chrono>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <limits>
#include <vector>

using namespace std;

void
sortit(vector<unsigned> &x)
{
    sort(x.begin(), x.end());
}

vector<unsigned>
unique_fill(gsl_rng const *r, unsigned n)
{
    vector<unsigned> x;
    for (unsigned i = 0; i < n; ++i)
        {
            auto t = static_cast<unsigned>(gsl_ran_flat(r, 0, 10000));
            while (find(x.begin(), x.end(), t) < x.end())
                {
                    t = static_cast<unsigned>(gsl_ran_flat(r, 0, 10000));
                }
            x.push_back(t);
        }
    sortit(x);
    return x;
}

void
mutrec_current(vector<unsigned> &recombinant, const vector<unsigned> &g1,
               const vector<unsigned> &g2, const vector<unsigned> &muts,
               const vector<unsigned> &brk)
{
    // Do recombination a la fwdpp
    // vector<unsigned> recombinant;
    recombinant.clear();
    auto b1 = g1.begin();
    auto e1 = g1.end();
    auto b2 = g2.begin();
    auto e2 = g2.end();

    for (auto &&b : brk)
        {
            auto itr = upper_bound(b1, e1, b);
            recombinant.insert(recombinant.end(), b1, itr);
            b1 = itr;
            b2 = upper_bound(b2, e2, b);
            swap(b1, b2);
            swap(e1, e2);
        }
    // Add mutation: this is (relatively)
    // inefficient
    for (auto &&m : muts)
        {
            recombinant.insert(
                upper_bound(recombinant.begin(), recombinant.end(), m), m);
        }
}

template <typename itr>
itr
update_mutation_iterator(std::vector<unsigned> &recombinant, itr mb, itr me,
                         itr b1, itr e1, unsigned b)
{
    while (mb < me && b1 < e1 && *mb < b && *mb < *b1)
        {
            recombinant.push_back(*mb);
            ++mb;
        }
    return mb;
}

void
mutrec_new(vector<unsigned> &recombinant, const vector<unsigned> &g1,
           const vector<unsigned> &g2, const vector<unsigned> &muts,
           const vector<unsigned> &brk)
{
    recombinant.clear();
    auto b1 = g1.begin();
    auto e1 = g1.end();
    auto b2 = g2.begin();
    auto e2 = g2.end();
    auto mb = muts.begin();
    auto me = muts.end();

    // for (auto b = brk.begin(); b != brk.end();
    //     ++b) // iterate over crossover positions
    for (auto &&b : brk)
        {
            // Take any mutations with positions < *b && < *b1
            // and add them to the recombinant
            mb = update_mutation_iterator(recombinant, mb, me, b1, e1, b);
            while (mb < me && *mb < b)
                {
                    auto itr = upper_bound(b1, e1, *mb);
                    recombinant.insert(recombinant.end(), b1, itr);
                    recombinant.push_back(*mb);
                    b1 = itr;
                    // mb = update_mutation_iterator(recombinant, ++mb, me, b1,
                    // e1, b);
                    ++mb;
                }
            auto itr = upper_bound(b1, e1, b);
            recombinant.insert(recombinant.end(), b1, itr);
            b1 = itr;
            b2 = upper_bound(b2, e2, b);
            swap(b1, b2);
            swap(e1, e2);
        }
    recombinant.insert(recombinant.end(), mb, me);
}

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
    chrono::duration<double> td1(0.), td2(0.);
    vector<unsigned> recombinant1,recombinant2;
    recombinant1.reserve((mg1 + mg2) / 2);
    recombinant2.reserve((mg1 + mg2) / 2);
    for (unsigned i = 0; i < 50000; ++i)
        {
            auto g1 = unique_fill(r, gsl_ran_poisson(r, mg1));
            auto g2 = unique_fill(r, gsl_ran_poisson(r, mg2));
            auto muts = unique_fill(r, gsl_ran_poisson(r, mnm));
            auto brk = unique_fill(r, gsl_ran_poisson(r, mbrk));
            brk.push_back(numeric_limits<unsigned>::max());
            std::chrono::time_point<std::chrono::system_clock> start, end,
                start2, end2;
            start = std::chrono::system_clock::now();
            mutrec_new(recombinant1, g1, g2, muts, brk);
            end = std::chrono::system_clock::now();
            start2 = std::chrono::system_clock::now();
            mutrec_current(recombinant2, g1, g2, muts, brk);
            end2 = std::chrono::system_clock::now();
            chrono::duration<double> d1 = end - start, d2 = end2 - start2;
            td1 += d1;
            td2 += d2;
            if (recombinant1 != recombinant2)
                {
                    cerr << "failure after " << i << " successes\n";
                }
        }
    cout << td1.count() / 50000. << ' ' << td2.count() / 50000. << '\n';
}
