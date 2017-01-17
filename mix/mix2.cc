#include <algorithm>
#include <chrono>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <limits>
#include <vector>
#include <queue>

using namespace std;

enum class NEXT_EVENT_TYPE
{
    MUT,
    REC,
    NONE
};

inline NEXT_EVENT_TYPE//std::pair<NEXT_EVENT_TYPE, double>
get_next_pos2(vector<unsigned>::iterator rb, vector<unsigned>::iterator re,
              vector<unsigned>::iterator mb, vector<unsigned>::iterator me)
{
    //if (rb == re && mb == me)
    //    {
	//		return NEXT_EVENT_TYPE::NONE;
    //        //return { NEXT_EVENT_TYPE::NONE,
    //        //         std::numeric_limits<unsigned>::quiet_NaN() };
    //    }
    if (rb == re || (mb!=me && *mb<*rb))
        {
			return NEXT_EVENT_TYPE::MUT;
            //return { NEXT_EVENT_TYPE::MUT, *mb };
        }
	else if (rb!=re) return NEXT_EVENT_TYPE::REC;

	return NEXT_EVENT_TYPE::NONE;
	/*
    if (mb == me)
        {
			return NEXT_EVENT_TYPE::REC;
            //return { NEXT_EVENT_TYPE::REC, *rb };
        }
    if (*mb < *rb)
        {
			return NEXT_EVENT_TYPE::MUT;
            //return { NEXT_EVENT_TYPE::MUT, *mb };
        }
		*/
	//return NEXT_EVENT_TYPE::REC;
    //return { NEXT_EVENT_TYPE::REC, *rb };
}

inline std::pair<NEXT_EVENT_TYPE, double>
get_next_pos(const std::queue<unsigned> &rec, const std::queue<unsigned> &mut)
{
    if (rec.empty() && mut.empty())
        {
            return { NEXT_EVENT_TYPE::NONE,
                     std::numeric_limits<unsigned>::quiet_NaN() };
        }
    if (rec.empty())
        {
            return { NEXT_EVENT_TYPE::MUT, mut.front() };
        }
    if (mut.empty())
        {
            return { NEXT_EVENT_TYPE::REC, rec.front() };
        }
    if (mut.front() < rec.front())
        {
            return { NEXT_EVENT_TYPE::MUT, mut.front() };
        }
    return { NEXT_EVENT_TYPE::REC, rec.front()};
}

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

void
mutrec_new2(vector<unsigned> &recombinant, const vector<unsigned> &g1,
            const vector<unsigned> &g2, vector<unsigned> &muts,
            vector<unsigned> &brk)
{
    recombinant.clear();
    auto b1 = g1.begin();
    auto e1 = g1.end();
    auto b2 = g2.begin();
    auto e2 = g2.end();
    auto mb = muts.begin(), me = muts.end(), rb = brk.begin(), re = brk.end();
    // while (!muts.empty() || !brk.empty())
    while (rb < re || mb < me)
        {
            auto event = get_next_pos2(rb, re, mb, me);
            if (event == NEXT_EVENT_TYPE::MUT)
                {
                    auto itr = upper_bound(b1, e1, *mb);
                    recombinant.insert(recombinant.end(), b1, itr);
                    recombinant.push_back(*mb);
                    b2 = upper_bound(b2, e2, *mb);
                    b1 = itr;
					++mb;
                }
            else if (event== NEXT_EVENT_TYPE::REC)
                {
                    auto itr = upper_bound(b1, e1, *rb);
                    recombinant.insert(recombinant.end(), b1, itr);
                    b2 = upper_bound(b2, e2, *rb);
                    b1 = itr;
                    swap(b1, b2);
                    swap(e1, e2);
					++rb;
                }
            else
                {
                    recombinant.insert(recombinant.end(), b1, e1);
                }
        }
}

void
mutrec_new(vector<unsigned> &recombinant, const vector<unsigned> &g1,
           const vector<unsigned> &g2, queue<unsigned> &muts,
           queue<unsigned> &brk)
{
    recombinant.clear();
    auto b1 = g1.begin();
    auto e1 = g1.end();
    auto b2 = g2.begin();
    auto e2 = g2.end();
    while (!muts.empty() || !brk.empty())
        {
            // std::cerr << std::distance(mb,me) << ' ' << std::distance(rb,re)
            // << '\n';
            auto event = get_next_pos(brk, muts);
            if (event.first == NEXT_EVENT_TYPE::MUT)
                {
                    auto itr = upper_bound(b1, e1, event.second);
                    recombinant.insert(recombinant.end(), b1, itr);
                    recombinant.push_back(event.second);
                    b2 = upper_bound(b2, e2, event.second);
                    b1 = itr;
					muts.pop();
                }
            else if (event.first == NEXT_EVENT_TYPE::REC)
                {
                    auto itr = upper_bound(b1, e1, event.second);
                    recombinant.insert(recombinant.end(), b1, itr);
                    b2 = upper_bound(b2, e2, event.second);
                    b1 = itr;
                    swap(b1, b2);
                    swap(e1, e2);
					brk.pop();
                }
            else
                {
                    recombinant.insert(recombinant.end(), b1, e1);
                }
        }
    recombinant.insert(recombinant.end(), b1, e1);
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
    vector<unsigned> recombinant1, recombinant2;
    recombinant1.reserve((mg1 + mg2) / 2);
    recombinant2.reserve((mg1 + mg2) / 2);
    unsigned successes = 0;
    for (unsigned i = 0; i < 50000; ++i)
        {
            auto g1 = unique_fill(r, gsl_ran_poisson(r, mg1));
            auto g2 = unique_fill(r, gsl_ran_poisson(r, mg2));
            auto muts = unique_fill(r, gsl_ran_poisson(r, mnm));
            auto brk = unique_fill(r, gsl_ran_poisson(r, mbrk));
            queue<unsigned> muts2(
                std::deque<unsigned>(muts.begin(), muts.end())),
                brk2(std::deque<unsigned>(brk.begin(), brk.end()));
            brk.push_back(numeric_limits<unsigned>::max());
            std::chrono::time_point<std::chrono::system_clock> start, end,
                start2, end2;
            start = std::chrono::system_clock::now();
            //mutrec_new(recombinant1, g1, g2, muts2, brk2);
            mutrec_new2(recombinant1, g1, g2, muts, brk);
            end = std::chrono::system_clock::now();
            start2 = std::chrono::system_clock::now();
            mutrec_current(recombinant2, g1, g2, muts, brk);
            end2 = std::chrono::system_clock::now();
            chrono::duration<double> d1 = end - start, d2 = end2 - start2;
            td1 += d1;
            td2 += d2;
            if (recombinant1 != recombinant2)
                {
                    cerr << "failure after " << successes << " successes and "
                         << i << " attempts: " << recombinant1.size() << ' '
                         << recombinant2.size() << '\n';
                    for (auto &&gi : g1)
                        cerr << gi << ' ';
                    cerr << '\n';
                    for (auto &&gi : g2)
                        cerr << gi << ' ';
                    cerr << "\nmuts: ";
                    for (auto &&mi : muts)
                        std::cerr << mi << ' ';
                    cerr << "\nbrk: ";
                    for (auto &&bi : brk)
                        cerr << bi << ' ';
                    std::cerr << "\n(";
                    for (auto &&ri : recombinant1)
                        std::cerr << ri << ' ';
                    std::cerr << ")\n(";
                    for (auto &&ri : recombinant2)
                        std::cerr << ri << ' ';
                    std::cerr << ")\n";
                }
            else
                {
                    ++successes;
                }
        }
    cout << td1.count() / 50000. << ' ' << td2.count() / 50000. << '\n';
}
