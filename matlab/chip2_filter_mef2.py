import sys

#l = sys.stdin.next()
#print l,

chip_ind = range(1, 16, 3)
q_ind = range(3, 16, 3)

for l in sys.stdin:
    t = l.strip().split('\t')
    try:
	qvals = map(float, map(t.__getitem__, q_ind))
	chipvals = map(float, map(t.__getitem__, chip_ind))
	qcond = map(lambda x: x < 1.0, qvals)
	chipcond = map(lambda x: x > 0.7, chipvals)
	fullcond = map(lambda x, y: x and y, qcond, chipcond)
	if reduce(lambda x, y: x or y, fullcond, False):
	    #print '\t'.join(t)
	    print '%s\t%s' % (t[0], '\t'.join(map(str, map(int, fullcond))))
    except:
	pass
