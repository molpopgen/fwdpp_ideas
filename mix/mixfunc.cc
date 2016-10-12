#include <vector>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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
mutrec_current(const vector<unsigned> &g1, const vector<unsigned> &g2,
               const vector<unsigned> &muts, const vector<unsigned> &brk)
{
    // for (auto &&i : g1)
    // cout << i << ' ';
    // cout << '\n';
    // for (auto &&i : g2)
    // cout << i << ' ';
    // cout << '\n';
    // for (auto &&i : muts)
    // cout << i << ' ';
    // cout << '\n';
    // for (auto &&i : brk)
    // cout << i << ' ';
    // cout << '\n';

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
    // for (auto &&i : recombinant)
    // cout << i << ' ';
    // cout << '\n';
    return recombinant;
}

template <typename itr>
itr
update_mutation_iterator(std::vector<unsigned> &recombinant, itr mb, itr me,
                         itr b1, itr e1, itr b)
{
    while (mb < me && b1 < e1 && *mb < *b && *mb < *b1)
        {
            recombinant.push_back(*mb);
            ++mb;
        }
    return mb;
}
vector<unsigned>
mutrec_new(const vector<unsigned> &g1, const vector<unsigned> &g2,
           const vector<unsigned> &muts, const vector<unsigned> &brk)
{
    // cout << "g1: ";
    // for (auto &&i : g1)
    // cout << i << ' ';
    // cout << '\n';
    // cout << "g2: ";
    // for (auto &&i : g2)
    // cout << i << ' ';
    // cout << '\n';
    // cout << "muts: ";
    // for (auto &&i : muts)
    // cout << i << ' ';
    // cout << '\n';
    // cout << "breakpoints: ";
    // for (auto &&i : brk)
    // cout << i << ' ';
    // cout << '\n';
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
            mb = update_mutation_iterator(recombinant, mb, me, b1, e1, b);
            while (mb < me && *mb < *b)
                {
                    auto itr = upper_bound(b1, e1, *mb);
                    recombinant.insert(recombinant.end(), b1, itr);
                    recombinant.push_back(*mb);
                    b1 = itr;
					//mb = update_mutation_iterator(recombinant, ++mb, me, b1, e1, b);
                    ++mb;
                }
            auto itr = upper_bound(b1, e1, *b);
            recombinant.insert(recombinant.end(), b1, itr);
            b1 = itr;
            b2 = upper_bound(b2, e2, *b);
            swap(b1, b2);
            swap(e1, e2);
        }
    recombinant.insert(recombinant.end(), mb, me);
    // cout << "result: ";
    // for (auto &&i : recombinant)
    // cout << i << ' ';
    // cout << '\n';
    return recombinant;
}
