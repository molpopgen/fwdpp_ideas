#Random ideas to try

## Object recycling

Memory operations are still a large fraction of run time.  Instead of removing extinct stuff, we can instead add a flag for "recycled" to mutations and gametes.  When we get a new object, we can simply update the first "recycled" object in the container, if one exists, else we can insert a new one.