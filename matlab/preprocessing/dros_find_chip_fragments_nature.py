import sys, os, numpy

DATADIR = './'

def binary_search(a, x, lo=0, hi=None):
    if hi is None:
        hi = len(a)
    while lo < hi:
        mid = (lo+hi)//2
        if x < a[mid][0]: hi = mid
        else: lo = mid+1
    return lo

class genome_dataset:
    def matches_from_point(self, chromosome, point, distance):
        startind = binary_search(self.fragments[chromosome], point-distance)
        stopind = binary_search(self.fragments[chromosome], point+distance)
        return self.fragments[chromosome][startind:stopind]

    def matches_from_gene(self, gene, distance):
        [chromosome, start, stop] = genes.genes[gene][0:3]
        #print chromosome, start, stop
        startind = binary_search(self.fragments[chromosome], start-distance)
        stopind = binary_search(self.fragments[chromosome], stop+distance)
        #print startind, stopind
        return self.fragments[chromosome][startind:stopind]

    def distance_to_closest_match(self, gene):
        [chromosome, start, stop] = genes.genes[gene][0:3]
        #print chromosome, start, stop
        startind = binary_search(self.fragments[chromosome], start)
        # Gene is before any fragment
        if startind == 0:
            frag = self.fragments[chromosome][startind]
            return max(frag[0] - stop, 0)
        # Usual case: gene in the middle of fragments
        elif startind < len(self.fragments[chromosome]):
            [frag1, frag2] = self.fragments[chromosome][(startind-1):(startind+1)]
            #print frag1, frag2
            return max(min(start - frag1[1], frag2[0] - stop), 0)
        # Gene is after all the fragments
        else:
            frag = self.fragments[chromosome][-1]
            return max(start - frag[1], 0)
        
class dros_genes(genome_dataset):
    def __init__(self):
        f = open(DATADIR + 'dm_genes.txt')
        f.next()
        self.fragments = dict()
        self.genes = dict()
        for l in f:
            t = l.strip().split('\t')
            if t[1] not in self.fragments:
                self.fragments[t[1]] = [map(int, t[2:4]) + [t[0]]]
            else:
                self.fragments[t[1]].append(map(int, t[2:4]) + [t[0]])
            self.genes[t[0]] = [t[1]] + map(int, t[2:4])
        for k in self.fragments.values():
            k.sort()

genes = dros_genes()

class nature_fragments(genome_dataset):
    def __init__(self, tf):
        f = open(DATADIR + 'nature08531-s3.txt')
        for i in range(6):
            l = f.next()
        header = l.strip().split('\t')
        indices = numpy.array(map(lambda x: tf in x, header))
        self.fragments = dict()
        for l in f:
            t = numpy.array(l.strip().split('\t'))
            if any(map(float, t[indices])):
                if t[1] not in self.fragments:
                    self.fragments[t[1]] = [map(float, t[2:4]) + [t[0]]]
                else:
                    self.fragments[t[1]].append(map(float, t[2:4]) + [t[0]])
        for k in self.fragments.values():
            k.sort()



tfs = ['mef2', 'twi']

frags = []
for k in range(len(tfs)):
    frags.append(dros_genome_datasets.nature_fragments(tfs[k]))

genes = dros_genome_datasets.genes

distances = dict()

for gene in genes.genes.keys():
    #print gene
    distances[gene] = map(lambda x: x.distance_to_closest_match(gene), frags)

v = distances.items()
v.sort()
print 'FBgn\t' + '\t'.join(map(lambda x: x+'_distance', tfs))
for t in v:
    print '%s\t%d\t%d' % ((t[0], ) + tuple(map(int, t[1])))
