#!/usr/bin/env python2.7

import unittest
from dnaseqlib import *

### Utility classes ###

# Maps integer keys to a set of arbitrary values.
class Multidict:
    # Initializes a new multi-value dictionary, and adds any key-value
    # 2-tuples in the iterable sequence pairs to the data structure.
    def __init__(self, pairs=[]):
        self.multidict = dict()
        for k, v in pairs:
            self.put(k, v)
    # Associates the value v with the key k.
    def put(self, k, v):
        if k in self.multidict:
            self.multidict[k].append(v)
        else:
            self.multidict[k] = [v]
    # Gets any values that have been associated with the key k; or, if
    # none have been, returns an empty sequence.
    def get(self, k):
       return self.multidict.get(k, [])


# Given a sequence of nucleotides, return all k-length subsequences
# and their hashes.  (What else do you need to know about each
# subsequence?)
def subsequenceHashes(seq, k):
    seq = iter(seq)  

    subseq = ''
    for i in range(k):
        c = next(seq, None)
        if c is None:
            return  
        subseq += c

    rh = RollingHash(subseq)
    pos = 0
    yield (rh.current_hash(), (pos, subseq))

    while True:
        c = next(seq, None)
        if c is None:
            break
        prev = subseq[0]
        subseq = subseq[1:] + c
        rh.slide(prev, c)
        pos += 1
        yield (rh.current_hash(), (pos, subseq))



# Similar to subsequenceHashes(), but returns one k-length subsequence
# every m nucleotides.  (This will be useful when you try to use two
# whole data files.)
def intervalSubsequenceHashes(seq, k, m):
    seq = iter(seq) 

    # Get the first k-length substring to initialize
    subseq = ''
    for _ in range(k):
        c = next(seq, None)
        if c is None:
            return  
        subseq += c

    rh = RollingHash(subseq) 
    pos = 0
    yield (rh.current_hash(), (pos, subseq))  

    # Then slide forward one char at a time, but only yield every m steps
    while True:
        try:
            for _ in range(m):  
                prev = subseq[0] 
                c = next(seq)    
                subseq = subseq[1:] + c
                rh.slide(prev, c)
                pos += 1
            yield (rh.current_hash(), (pos, subseq)) 
        except StopIteration:
            return  # Stop when sequence is exhausted
       

# Searches for commonalities between sequences a and b by comparing
# subsequences of length k.  The sequences a and b should be iterators
# that return nucleotides.  The table is built by computing one hash
# every m nucleotides (for m >= k).
def getExactSubmatches(a, b, k, m):
    table = Multidict()
    for h, (apos, asubseq) in intervalSubsequenceHashes(iter(a), k, m):
        table.put(h, (apos, asubseq))
    
    for h, (bpos, bsubseq) in subsequenceHashes(iter(b), k):
        for apos, asubseq in table.get(h):
            if asubseq == bsubseq:
                yield (apos, bpos)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Usage: {0} [file_a.fa] [file_b.fa] [output.png]'.format(sys.argv[0]))
        sys.exit(1)

    # The arguments are, in order: 1) Your getExactSubmatches
    # function, 2) the filename to which the image should be written,
    # 3) a tuple giving the width and height of the image, 4) the
    # filename of sequence A, 5) the filename of sequence B, 6) k, the
    # subsequence size, and 7) m, the sampling interval for sequence
    # A.
    compareSequences(getExactSubmatches, sys.argv[3], (500,500), sys.argv[1], sys.argv[2], 8, 100)
