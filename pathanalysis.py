#!/usr/bin/env python
# coding: utf-8
#searchの出力を読みこみ、始点(原点)から最短距離xにある点に行く、他の経路の本数の統計をとる。
#各行の最後の格子点をまず読む。
#ループがないことを確認する(同じ数字がでてこない=SAWである)
#最後の格子点の統計として記録する。

#import sys


sites = []
sitefile = open("60x1234.ar3a")
while True:
    line = sitefile.readline()
    if line == "":
        break
    if line[0:5] == "@AR3A":
        line = sitefile.readline()
        nsite = int(line)
        for i in range(nsite):
            line = sitefile.readline()
            columns = line.split()
            sites.append(map(float,columns))

def dictofdict_add(d,key,value):
    if not d.has_key(key):
        d[key] = dict()
    d[key][value] = 1

bondpos = dict()
bondowners = dict()
nodeowners = dict()
bondfile = open("60x1234.ngph")
while True:
    line = bondfile.readline()
    if line == "":
        break
    if line[0:5] == "@NGPH" or line[0:5] == "@PPOS":
        line = bondfile.readline()
        nsite = int(line)
        while True:
            line = bondfile.readline()
            columns = line.split()
            columns = map(int,columns)
            if columns[0] < 0:
                break
            i,j,k = columns
            x = (sites[i][0] + sites[j][0])/2
            y = (sites[i][1] + sites[j][1])/2
            z = (sites[i][2] + sites[j][2])/2
            bondpos[k] = (x,y,z,i,j)
            dictofdict_add(bondowners,k,i)
            dictofdict_add(bondowners,k,j)
            dictofdict_add(nodeowners,i,k)
            dictofdict_add(nodeowners,j,k)

#distance along the undirected graph from 48773
distances = dict()
queue = []
# total number of nth neighbor along the undirected graph
nnei = [0]*20

def BreadthFirstSearch(distances,queue,depth):
    newq = []
    for bond in queue:
        if distances.has_key(bond):
            continue
        distances[bond] = depth
        nnei[depth] += 1
        ends = bondowners[bond].keys()
        #print bond, ends,depth
        neibonds = nodeowners[ends[0]].keys() + nodeowners[ends[1]].keys()
        newq += neibonds
    return newq

queue.append(48773)
for i in range(16):
    queue = BreadthFirstSearch(distances,queue,i)













paths = dict()
dist  = dict()
pathfile = open("60x14x1234.log")
sawpath = [0]*20
while True:
    line = pathfile.readline()
    if line == "":
        break
    columns = line.split()
    if len(columns) > 3:
        d = columns[2]
        columns = tuple(columns[3:])
        flag = dict()
        saw = True #self-avoiding walk
        for column in columns:
            if flag.has_key(column):
                saw = False
            flag[column] = 1
        dest = int(columns[-1])
        pathlen = len(columns)
        dist[dest] = d
        if saw:
            sawpath[len(columns)] += 1
            if not paths.has_key(dest):
                paths[dest] = [0] * 20
            paths[dest][pathlen] += 1
print sawpath

#from math import sqrt
#x0,y0,z0,i0,j0 = bondpos[48773]
#for dest in paths:
#    x,y,z,i,j = bondpos[int(dest)]
#    x -= x0
#    y -= y0
#    z -= z0
#    print dest,sqrt(x**2+y**2+z**2),x,y,z,i0,j0,i,j,
#    for k in range(20):
#        print paths[dest][k],
#    print


#Table 1: Real distance and available path distance
sum = dict()
for dest in paths:
    if distances.has_key(dest):
        d = distances[dest] # distance along the undirected graph
        for k in range(20):
            if paths[dest][k] > 0:
                break
        #k is the shortest path along the directed graph
        if not sum.has_key(d):
            sum[d] = [0]*20
        sum[d][k] += 1
for d in sum:
    print d,
    wsum=0.0
    wei=0.0
    for dist in range(20):
        print sum[d][dist],
        wsum += sum[d][dist]*(dist-1)
        wei  += sum[d][dist]
    print wsum/wei,wei,nnei[d] #some paths are not taken into consideration when wei differs from nnei[d].

print

#Table 2: Number of detours
for dist in range(2,20):
    sum = [0.0] * 20
    cnt = [0.0] * 20
    for dest in paths:
        zero=True
        for k in range(dist):
            if paths[dest][k] > 0:
                zero = False
        if zero and paths[dest][dist]>0:
            #print dest,paths[dest]
            for k in range(20):
                sum[k] += paths[dest][k]
                cnt[k] += 1
    for k in range(20):
        if cnt[k]:
            print sum[k]/cnt[k],
        else:
            print "-",
    print

       
