# -*- coding: utf-8 -*-
"""
Created on September 1st 2015

@author: Nicolas Chrysanthos
"""

import numpy as np
import random
import math
import csv

rand_seed = 5674

n_seqs = 20 #number of sequences
lmin = 5 #minimum length of a sequence
lmax = 10 #maximum length of a sequence
ymin = 0 #min amplitude
ymax = 10 #max amplitude
sigma = 3. #bandwidth for the kernel

def createDummyData():
    """
    Creates a list of dummy sequences
    """

    random.seed(rand_seed)

    data = []

    for i in xrange(n_seqs):
        y0 = random.randrange(ymin, ymax)
        y1 = random.randrange(ymin, ymax)
        l = random.randint(lmin, lmax)

        seq = np.linspace(y0,y1, l)
        data.append(seq)
        #print seq

    return data

def reducedDistanceMatrix(seq1, seq2):
    """
    In the case of one-sided alignment, not all pairwise distances are needed,
    so we compute a reduced version of alignDistMat, of size n1*(n2 - n1 + 1)
    when n1<=n2.
    See Fig 4 in article.
    we return a matrix with vertical
    indexes always representing the *shorter* sequence
    """
    n1, n2 = len(seq1), len(seq2)
    
    if n1>n2:
        seqA, seqB = seq2, seq1
    else:
        seqA, seqB = seq1, seq2
    
    nA, nB = len(seqA), len(seqB)
    
    reducedDistMat = np.zeros((nA, nB - nA + 1))
    for i in xrange(nA):
        for j in xrange(nB - nA + 1):
            reducedDistMat[i,j] = (seqA[i] - seqB[j + i])**2 #replace with any distance that is Hilbertian
    
    return reducedDistMat#/nB #nB is the length of the *longer* sequence

def one_sided_mean(reducedDistMat):
    """
    Computes the mean of all One-Sided distances,
    with reducedDistMat containing pairwise *distances*
    """
    L, M = np.shape(reducedDistMat)
    dtwM = np.zeros((L + 1, M + 1))
    
    for m in xrange(M + 1):
        dtwM[0, m] = 0
    for l in xrange(L + 1):
        dtwM[l, 0] = 0
        
    dtwM[1,1] = reducedDistMat[0,0]

    for l in xrange(1,L+1):
        _l = float(l)
        for m in xrange(1, M+1):
            _m = float(m)
            if not (_l, _m)==(1,1):
                A = (_l - 1)/(_l + _m - 2) * dtwM[l-1, m]
                B = (_m - 1)/(_l + _m - 2) * dtwM[l, m-1]
                dtwM[l, m] = reducedDistMat[l - 1, m - 1] + A + B
    
    return dtwM[-1, -1]

if __name__ == '__main__':

    #Create dummy sequences
    #data = createDummyData()
    data = []
    with open('example_data.txt', 'rb') as csvfile:
        datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in datareader:
            row_float = [float(num) for num in row]
            data.append(row_float)

    #Now compute one by one the elements of the Gram matrix
    gramMat = np.zeros((n_seqs, n_seqs))
    for i in xrange(n_seqs):
        seq1 = data[i]
        l1 = len(seq1)
        for j in xrange(i,n_seqs):
            seq2 = data[j]
            l2 = len(seq2)
            reducedMat = reducedDistanceMatrix(seq1, seq2)
            mean_dist = one_sided_mean(reducedMat)/max(l1,l2) #divide by the length of the longer sequence
            gramMat[i,j] = math.exp(-mean_dist/(2*sigma)) #now just take the exponent with a bit of rescaling
            gramMat[j,i] = gramMat[i,j]

    #Print our beautiful matrix
    print "Gram Matrix:"
    print gramMat

    #Check that the matrix is psd
    d, v = np.linalg.eigh(gramMat)
    print "Eigenvalues:"
    print d