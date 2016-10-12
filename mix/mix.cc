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

vector<unsigned>
mutrec_new(gsl_rng const *r, const vector<unsigned> &g1,
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
    auto mb = muts.begin();
    auto me = muts.end();
    for (auto &&b : brk)
        {
            cout << "processing breakpoint : " << b << ":\n";
            auto itr = upper_bound(b1, e1, b);
            cout << "mutations in gamete to be inserted before " << b
                 << " are: ";
            for_each(b1, itr, [](unsigned i) { cout << i << ' '; });
            cout << '\n';
            cout << "new mutations with positions <= " << b << " are: ";
            auto mitr = upper_bound(mb, me, b);
            for_each(mb, mitr, [](unsigned i) { cout << i << ' '; });
            cout << '\n';
            cout << "new mutations that come before the gamete's mutations: "
                    "are: ";
            auto mb2 = mb;
            for (; mb2 < me && *mb2 < *b1 && *mb2 < b; ++mb2)
                {
                    cout << *mb2 << ' ';
                    recombinant.push_back(*mb2);
                }
            cout << '\n';
            cout << "new mutations interdigitated with the gamete's mutations "
                    "are: ";
            mitr = upper_bound(b1, itr, *mb2);
            while (mitr != b1 && mitr != itr)
                {
                    cout << *mb2 << ' ';
                    auto itr2 = upper_bound(b1, itr, *mb2);
                    cout << "(" << (itr2 == itr) <<','<<(itr2==b1)<< ','<<(b1<itr2)<<") ";
                    if (itr2 == b1)
                        {
                            recombinant.push_back(*mb2);
                        }
                    else
                        {
							for(;b1<itr2;++b1)
							{
								cout << '['<<*b1<<','<<*mb2<<']' << ' ';
								if(*b1<*mb2) recombinant.push_back(*b1);
							}
								recombinant.push_back(*mb2);
						}
                    ++mb2;
                    mitr = upper_bound(b1, itr, *mb2);
                }
            cout << '\n';
            //mb = upper_bound(mb2, me, b);
            itr = upper_bound(b1, e1, b);
            recombinant.insert(recombinant.end(), b1, itr);
            cout << "new mutations after gamete but before " << b << ": ";
			mb=mb2;
			for(;mb<me&&*mb<b;++mb)
			{
				cout << *mb << ' ';
				recombinant.push_back(*mb);
			}
            cout << '\n';
            b1 = itr;
            b2 = upper_bound(b2, e2, b);
            swap(b1, b2);
            swap(e1, e2);
        }
    for (auto &&i : recombinant)
        {
            if (binary_search(muts.begin(), muts.end(), i))
                {
                    cout << i << "* ";
                }
            else
                {
                    cout << i << ' ';
                }
        }
    cout << '\n';
    return recombinant;
}

int
main(int argc, char **argv)
{
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, 0);
    unsigned x = 5;
    for (unsigned i = 0; i < 3; ++i)
        {
            auto g1 = unique_fill(r, x);
            auto g2 = unique_fill(r, x);
            auto muts = unique_fill(r, x);
            auto brk = unique_fill(r, x);
            brk.push_back(numeric_limits<unsigned>::max());
            auto output1 = mutrec_current(r, g1, g2, muts, brk);
            auto output2 = mutrec_new(r, g1, g2, muts, brk);
            cout << "check: " << (output1 == output2) << '\n';
        }
}
