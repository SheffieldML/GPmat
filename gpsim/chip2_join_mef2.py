import sys, numpy

ind = {}
res = {}

f = open(sys.argv[1])
for l in f:
    t = l.strip().split(',')
    ind[t[0]] = t[2]
f.close()

for l in sys.stdin:
    t = l.strip().split('\t')
    if ind.has_key(t[0]):
	if res.has_key(ind[t[0]]):
	    res[ind[t[0]]] += numpy.array(map(int, t[1:]))
	else:
	    res[ind[t[0]]] = numpy.array(map(int, t[1:]))

a = res.items()
a.sort()
for t in a:
    print '%s\t%s' % (t[0], '\t'.join(map(str, list(t[1]))))
