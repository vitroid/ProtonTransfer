#!/usr/bin/env python

import sys

dmax = int(sys.argv[1])

paths = dict()
for line in sys.stdin:
   columns = line.split()
   if len(columns) > 3:
       depth = dmax - int(columns[0])
       last = columns[-1]
       columns = columns[3:-1]
       if not last in columns:
           print depth,last
           if not paths.has_key(last):
               paths[last] = [0]*(dmax+1)
           paths[last][depth] += 1

for last in paths:
    print last,
    for i in range(dmax+1):
        print paths[last][i],
    print



       
