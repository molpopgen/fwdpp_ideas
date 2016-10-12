#ifndef MIX_H_
#define MIX_H_
#include <vector>
#include <gsl/gsl_rng.h>
void
sortit(std::vector<unsigned> &x);

std::vector<unsigned>
unique_fill(gsl_rng const *r, unsigned n);

std::vector<unsigned>
mutrec_current(const std::vector<unsigned> &g1, const std::vector<unsigned> &g2,
               const std::vector<unsigned> &muts, const std::vector<unsigned> &brk);

std::vector<unsigned>
mutrec_new(const std::vector<unsigned> &g1, const std::vector<unsigned> &g2,
           const std::vector<unsigned> &muts, const std::vector<unsigned> &brk);
#endif
