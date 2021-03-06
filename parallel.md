# Notes on parallelizing forward simulations.

## Calculation of mean fitness.

A typical fwdpp-based simulation will loop over diploids, set gamete counts to zero, and fill a vector of fitnesses.
For example:

```c++
for( const auto & dip : diploids )
{
    gametes[dip.first].n=0;
    gametes[dip.second].n=0;
    fitnesses.push_back(fitness_function(dip,gametes,mutations));
}
```

For many cases, the above is trivial to parallelize asynchronously. We need to break the computation up into ranges,
apply the fitness function within each range, and update the (pre-allocated) fitnesses vector.  Here are two thoughts on
how to do this:

```c++
struct frange
// Abstract out application of fitness function over a range.
// "Trick" is that gametes and mutations will be std::cref's
// wrapping the underlying types, which explains the foo.get()
// when passing to the fitness function ff.
{
    template <typename itr, typename gcont, typename mcont, typename wfunc>
    inline std::pair<std::vector<double>, double>
    operator()(const gcont &gametes, const mcont &mutations, itr beg,
               itr end, const wfunc &ff)
    {
        std::vector<double> w;
        w.reserve(std::distance(beg, end));
        double sum = 0.0;
        for (; beg < end; ++beg)
            {
                w.push_back(ff(*beg, gametes.get(), mutations.get()));
                sum += w.back();
            }
        return std::make_pair(std::move(w), sum);
    }
};
```

```c++
struct frange2
// This version directly fills in the fitness vector.
// The potential issue with this approach is contention
// due to false sharing at the boundaries of ranges
// updated by different threads.
{
    template <typename itr, typename gcont, typename mcont, typename wfunc>
    inline double
    operator()(const gcont &gametes, const mcont &mutations,
               std::reference_wrapper<std::vector<double>> wref,
               std::size_t i, itr beg, itr end, const wfunc &ff)
    {
        auto & w = wref.get();
        double sum = 0.0;
        for (; beg < end; ++beg, ++i)
            {
                w[i] = ff(*beg, gametes.get(), mutations.get());
                sum += w[i];
            }
        return sum;
    }
};
```

*A priori*, the second version should be faster.  It should also do an ok job avoiding false sharing.  However, limited
testing suggests no real differences in performance.  The latter implementation is also "unsafe", in that there's no
range-checking when assiging into `w`.

I should read more about how to properly use the `std::reference_wrapper` types...

## Reproducibility 

A major challenge is how to get mutation generation to be reproducible AND parallelized.  For infinitely-many sites
models, the challenge is at the level of the lookup table.  A naive approach would be to protect access with some form
of spin lock and a relaxed memory order.  The result would be that simulation results are random for a given initial
seed, even if each thread had its on rng + seed.  The reason is that entries into the lookup table would be essentially
random, affecting the range of possible new mutation positions generated by subsequent threads.

It would seem that the solution involves more complex process management. If one process generated all mutations for all
descendants, the remaining processes could be generating recombination breakpoints for each new descendant.
