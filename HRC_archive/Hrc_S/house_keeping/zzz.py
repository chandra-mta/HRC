#!/usr/bin/env /proj/sot/ska/bin/python

f = open('cancelled_list', 'r')
data = [line.strip() for line in f.readlines()]
f.close()

sdata = list(set(data))

fo = open('xxx', 'w')
for ent in sdata:
    fo.write(ent)
    fo.write('\n')
fo.close()
