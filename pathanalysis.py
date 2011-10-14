#!/usr/bin/env python
# coding: utf-8

#work/*.shの出力から、再帰確率を算出する。


import sys


p = [[0.0]*20 for i in range(20)]
n = [0.0] * 20
for line in sys.stdin:
    columns = line.split()
    columns = map(float, columns[1:])
    for i in range(len(columns)):
        if columns[i]> 0:
            break
    n[i] += 1
    for j in range(i,len(columns)):
        p[i][j] += columns[j]


for i in range(20):
    if n[i]:
        #print i,
        for j in range(20):
            p[i][j]/=n[i]
        #    print p[i][j],
        #print

for i in range(20):
    d = 1.0
    sum = 0.0
    for j in range(20):
       #print p[i][j]/d,
       sum += p[i][j]/d
       print i,j,sum
       d = 4*(d-p[i][j])
       if d == 0.0:
           break
    print

