#!/usr/bin/env python
from Factor import *
import numpy as np
from itertools import product
import sys
import itertools

def getUniqueVar( factorList):
    """ given factorList which is a list of Factor objects
        return a list of unique variables appearing in the
        Factor objects. See http://stackoverflow.com/a/2151553/1735942 """

    return list(set().union(*[  list(f.getVar())  for f in factorList ] ))


def isMemberBoolean (A, B):
    """  returns an list of the same  as A containing True where the elements of A are in B and False otherwise """
    
    return [ x in B for x in A ]

def isMember( A, B):
    """ return a python list containing  indices in B where the elements of A are located
        A and B are numpy 1-d arrays
        mapA[i]=j if and only if B[i] == A[j]"""
    mapA=[]
    for i in range(len(A)):
        mapA.append( np.where(B==A[i])[0].tolist()[0] )

    return mapA




def accum(accmap, a, func=None, size=None, fill_value=0, dtype=None):
    """ based on the recipie here to implement a numpy equivalent to matlab accumarray:
        http://www.scipy.org/Cookbook/AccumarrayLike
        An accumulation function similar to Matlab's `accumarray` function.

    Parameters
    ----------
    accmap : ndarray
        This is the "accumulation map".  It maps input (i.e. indices into
        `a`) to their destination in the output array.  The first `a.ndim`
        dimensions of `accmap` must be the same as `a.shape`.  That is,
        `accmap.shape[:a.ndim]` must equal `a.shape`.  For example, if `a`
        has shape (15,4), then `accmap.shape[:2]` must equal (15,4).  In this
        case `accmap[i,j]` gives the index into the output array where
        element (i,j) of `a` is to be accumulated.  If the output is, say,
        a 2D, then `accmap` must have shape (15,4,2).  The value in the
        last dimension give indices into the output array. If the output is
        1D, then the shape of `accmap` can be either (15,4) or (15,4,1)
    a : ndarray
        The input data to be accumulated.
    func : callable or None
        The accumulation function.  The function will be passed a list
        of values from `a` to be accumulated.
        If None, numpy.sum is assumed.
    size : ndarray or None
        The size of the output array.  If None, the size will be determined
        from `accmap`.
    fill_value : scalar
        The default value for elements of the output array.
    dtype : numpy data type, or None
        The data type of the output array.  If None, the data type of
        `a` is used.

    Returns
    -------
    out : ndarray
        The accumulated results.

        The shape of `out` is `size` if `size` is given.  Otherwise the
        shape is determined by the (lexicographically) largest indices of
        the output found in `accmap`.


    Examples
    --------
     from numpy import array, prod
     a = array([[1,2,3],[4,-1,6],[-1,8,9]])
     a
    array([[ 1,  2,  3],
           [ 4, -1,  6],
           [-1,  8,  9]])
     # Sum the diagonals.
     accmap = array([[0,1,2],[2,0,1],[1,2,0]])
     s = accum(accmap, a)
    array([9, 7, 15])
     # A 2D output, from sub-arrays with shapes and positions like this:
     # [ (2,2) (2,1)]
     # [ (1,2) (1,1)]
     accmap = array([
            [[0,0],[0,0],[0,1]],
            [[0,0],[0,0],[0,1]],
            [[1,0],[1,0],[1,1]],
        ])
     # Accumulate using a product.
     accum(accmap, a, func=prod, dtype=float)
    array([[ -8.,  18.],
           [ -8.,   9.]])
     # Same accmap, but create an array of lists of values.
      accum(accmap, a, func=lambda x: x, dtype='O')
      array([[[1, 2, 4, -1], [3, 6]],
           [[-1, 8], [9]]], dtype=object)
    """

    if accmap.shape[:a.ndim] != a.shape:
        raise ValueError("The initial dimensions of accmap must be the same as a.shape")
    if func is None:
        func = np.sum
    if dtype is None:
        dtype = a.dtype
    if accmap.shape == a.shape:
        accmap = np.expand_dims(accmap, -1)
    adims = tuple(range(a.ndim))
    if size is None:
        size = 1 + np.squeeze(np.apply_over_axes(np.max, accmap, axes=adims))
        size = np.atleast_1d(size)
    size=np.array(size, dtype=int)
    # Create an array of python lists of values.
    vals = np.empty(size, dtype='O')
    for s in product(*[range(k) for k in size]):
        vals[s] = []
    for s in product(*[range(k) for k in a.shape]):
        indx = tuple(accmap[s])
        val = a[s]
        vals[indx].append(val)

   # Create the output array.
    out = np.empty(size, dtype=dtype)
    for s in product(*[range(k) for k in size]):
        if vals[s] == []:
            out[s] = fill_value
        else:
            out[s] = func(vals[s])

    return out



def getIndex ( V, I):
    """ this method finds for every element of 1d- NumPy array, the index in another 1d-NumPy array
        based off of this StackOverflow answer http://stackoverflow.com/a/8251757/1735942
        V is a numpy 1d array listing the variables of a factor in a PGM
        I is the numpy list reprsenting the intersection of variables between two factors

        This method returns the index of the intersection variables in the list V 
        """
    index=np.argsort(V) #np.argsort gives the index ordering of an array that would result in sorted order
                        #http://docs.scipy.org/doc/numpy/reference/generated/numpy.argsort.html
    sorted_v=V[index]   #sorted version of V array
    sorted_index=np.searchsorted( sorted_v, I) #i ndices where elements of I  would be inserted into sorted_v to maintain order
                                               #http://docs.scipy.org/doc/numpy/reference/generated/numpy.searchsorted.html
    vindex=np.take(index, sorted_index, mode='clip') #Take elements from an array along an axis.\
                                                     #http://docs.scipy.org/doc/numpy/reference/generated/numpy.take.html
    #print vindex
    return vindex.tolist()

def generateAlleleGenotypeMappers( numAlleles):
    """ analogous to AssignmentToIndex and IndexToAssignment
        this function maps alleles to genotypes and genotypes to alleles
        each allele's id is its index in the list of alleles/list of allele freqs
        Similar for genotypes, its id is its index in list of genotypes

        This function returns to numpy 2-d arrays. if n is the number of alleles
        and m is the number of genotypes, then the following NumPy data structures
        are returned

        allelesToGenotypes: n x n matrix that maps pairs of allele IDs to
        genotype IDs -- if allelesToGenotypes[i, j] = k, then the genotype with
        ID k comprises of the alleles with IDs i and j

        genotypesToAlleles: m x 2 matrix of allele IDs, where m is the number of
        genotypes -- if genotypesToAlleles(k, :) = i, j, then the genotype with ID k
        is comprised of the allele with ID i and the allele with ID j

    """
    allelesToGenotypes=np.zeros([numAlleles,numAlleles], dtype=np.int)
    index=0
    for i in range(numAlleles):
        for j in range (i,numAlleles):
            allelesToGenotypes[i, j] = index;
            index+=1

    for i in range(numAlleles):
        for j in range(0,i):
            allelesToGenotypes[i,j] = allelesToGenotypes[j,i]


    numGenotypes= (numAlleles * (numAlleles - 1))/2 + numAlleles
    genotypesToAlleles=np.zeros([numGenotypes,2], dtype=np.int)

    index=0
    for i in range(numGenotypes):
        for j in range(i,numAlleles):
        #print i, j
            genotypesToAlleles[index, :] = [i, j];
            index+=1
    return ( allelesToGenotypes, genotypesToAlleles)


def genotypeToIndex( geno, alleles='ACGT'):
    """ given a string enumerating possible alleles
        return the index of geno in enumerated list of genotypes
        by default, the enumerated list is all 10 possibel genotypes"""

    genotypes= [ "".join(list(genotype))  for genotype in itertools.combinations_with_replacement(alleles, 2) ]
    print genotypes
    
    try:
        return  genotypes.index(geno)
    except ValueError:
        print "genotype not in list of genotypes."
    

def indexToGenotype( index, alleles='ACGT' ):
    """ return genotype at a given index position after
    enumerating all possible genotypes given string of alleles and
    assigning to a list. By default the list contains all possible 10 genotypes"""

    
    genotypes= [ "".join(list(genotype))  for genotype in itertools.combinations_with_replacement(alleles, 2) ]
    
    try:
       return genotypes[index]
    except IndexError:
        print "Index out of bounds, not a valid index for list of genotypes"