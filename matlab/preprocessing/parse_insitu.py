import re
import sys, os
import numpy

TERMPATH=os.path.expanduser('~/data/eileen_data/magnus-data/insitu')
TERMFILES = ['nature_terms_meso.txt', 'nature_terms_vm.txt',
             'nature_terms_sm.txt', 'nature_terms_cm.txt']

termsets = []

for k in range(len(TERMFILES)):
    f = open(os.path.join(TERMPATH, TERMFILES[k]))
    termsets.append(set())
    for l in f:
        #print '%d %s' % (k, l.strip())
        termsets[k].add(l.strip())
    f.close()

#print termsets

genestore = dict()
condstore = []

for k in range(len(termsets) + 1):
    condstore.append(dict())

# 1 - gene, 3 - annotation
sys.stdin.next()
for l in sys.stdin:
    t = l.strip().split('\t')
    if len(t[1])==0:
        continue
    #print t
    loc = []
    for k in range(len(termsets)):
        if t[3] in termsets[k]:
            #print '%s: %d' % (t[3], k)
            loc.append(k)

    if len(loc) == 0:
        if not genestore.has_key(t[1]):
            genestore[t[1]] = []
    for k in loc:
        if not genestore.has_key(t[1]):
            genestore[t[1]] = [k]
        elif k not in genestore[t[1]]:
            genestore[t[1]].append(k)

#print unknown_terms
#print genestore

f = open('insitu_data.txt', 'w')
t = map(lambda x: x.split('.')[0].split('_')[2], TERMFILES)
f.write('FBgn ' + ' '.join(t) + '\n')
l = genestore.items()
l.sort()
for t in l:
    ar = numpy.zeros((len(termsets)))
    for k in t[1]:
        ar[k] = 1
    f.write(('%s' + ' %d' * len(ar.flatten()) + '\n') %
            ((t[0], ) + tuple(ar.flatten())))
f.close()
