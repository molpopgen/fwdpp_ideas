#Random ideas to try

## Object recycling

Memory operations are still a large fraction of run time.  Instead of removing extinct stuff, we can instead add a flag for "recycled" to mutations and gametes.  When we get a new object, we can simply update the first "recycled" object in the container, if one exists, else we can insert a new one.

This ended up being _slower_.  The overhead of creating lookup tables for "recycled" gametes + managing a larger linked list more than offset any potential benefits.

## Hash keys for gametes

key = sum(hash of mutations).

With a default that is sum( has(mut. position + mut. count) )

Why: operator== can be re-written as:

if(a.key==b.key() { 
//compare mutations (control for hash collisions
}
return false

This should be way faster than the default operator==

Update: actually, should just has the memory address of the mutation.  Guaranteed unique.

Hashing addresses is likely not possible: http://stackoverflow.com/questions/14167455/is-it-possible-to-hash-pointers-in-portable-c03-code

### Results

Hashing based on positions of mutations worked, but was actually a ton of overhead.  results: slower sims.  I have the code in the branch dev_hash on my server.

## More detailed lookup structure

Instead of

~~~{cpp}
map< int, multimap<int,iterator > >
~~~

Let's try this:

~~~{cpp}
//Lookup[# neutral mutations][# selected mutations][first neutral mutation position][all gametes with same first selected mutation position
map< int, multimap<int, map< double, multimap<double, iterator> > > >
~~~

## Revisiting gamete_base::operator==

Is a less naive option possible?

* http://letsalgorithm.blogspot.com/2012/02/intersecting-two-sorted-integer-arrays.html
* http://articles.leetcode.com/2010/03/here-is-phone-screening-question-from.html, esp comment by "surrender"
