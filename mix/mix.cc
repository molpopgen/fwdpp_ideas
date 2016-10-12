// Attempting to merged mutation and recombination into
// a single step

#include <algorithm>
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

vector<unsigned>
mutrec_current(gsl_rng const *r, const vector<unsigned> &g1,
               const vector<unsigned> &g2, const vector<unsigned> &muts,
               const vector<unsigned> &brk)
{
    for (auto &&i : g1)
        cout << i << ' ';
    cout << '\n';
    for (auto &&i : g2)
        cout << i << ' ';
    cout << '\n';
    for (auto &&i : muts)
        cout << i << ' ';
    cout << '\n';
    for (auto &&i : brk)
        cout << i << ' ';
    cout << '\n';

    // Do recombination a la fwdpp
    vector<unsigned> recombinant;
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
    for (auto &&i : recombinant)
        cout << i << ' ';
    cout << '\n';
    return recombinant;
}

template <typename itr>
itr
update_mutation_iterator(std::vector<unsigned> &recombinant, itr mb, itr me,
                         itr b1, itr e1,itr b)
{
    while (mb < me && b1 < e1 && *mb < *b && *mb < *b1)
        {
            recombinant.push_back(*mb);
            ++mb;
        }
    return mb;
}
vector<unsigned>
mutrec_new(gsl_rng const *r, const vector<unsigned> &g1,
           const vector<unsigned> &g2, const vector<unsigned> &muts,
           const vector<unsigned> &brk)
{
    cout << "g1: ";
    for (auto &&i : g1)
        cout << i << ' ';
    cout << '\n';
    cout << "g2: ";
    for (auto &&i : g2)
        cout << i << ' ';
    cout << '\n';
    cout << "muts: ";
    for (auto &&i : muts)
        cout << i << ' ';
    cout << '\n';
    cout << "breakpoints: ";
    for (auto &&i : brk)
        cout << i << ' ';
    cout << '\n';
    vector<unsigned> recombinant;
    auto b1 = g1.begin();
    auto e1 = g1.end();
    auto b2 = g2.begin();
    auto e2 = g2.end();
    auto mb = muts.begin();
    auto me = muts.end();

    for (auto b = brk.begin(); b != brk.end();
         ++b) // iterate over crossover positions
        {
            // Take any mutations with positions < *b && < *b1
            // and add them to the recombinant
			mb=update_mutation_iterator(recombinant,mb,me,b1,e1,b);
			while(mb<me&& *mb<*b)
			{
				auto itr = upper_bound(b1, e1, *mb);
				recombinant.insert(recombinant.end(), b1, itr);
				recombinant.push_back(*mb);
				b1 = itr;
				++mb;
			}
            auto itr = upper_bound(b1, e1, *b);
            recombinant.insert(recombinant.end(), b1, itr);
            b1 = itr;
            b2 = upper_bound(b2, e2, *b);
            swap(b1, b2);
            swap(e1, e2);
        }
    recombinant.insert(recombinant.end(),mb,me);
	cout << "result: ";
    for (auto &&i : recombinant)
        cout << i << ' ';
    cout << '\n';
    return recombinant;
}

int
main(int argc, char **argv)
{
	int argn=1;
	unsigned seed =static_cast<unsigned>(atoi(argv[argn++])); 
    unsigned x = static_cast<unsigned>(atoi(argv[argn++]));

    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed);
    for (unsigned i = 0; i < 300; ++i)
        {
            auto g1 = unique_fill(r, gsl_ran_poisson(r,1+x));
            auto g2 = unique_fill(r, gsl_ran_poisson(r,1+x));
            auto muts = unique_fill(r, gsl_ran_poisson(r,1+2*x));
            auto brk = unique_fill(r, gsl_ran_poisson(r,1+x));
            brk.push_back(numeric_limits<unsigned>::max());
            auto output1 = mutrec_current(r, g1, g2, muts, brk);
            auto output2 = mutrec_new(r, g1, g2, muts, brk);
            cout << "check: " << (output1 == output2) << '\n';
        }
}
