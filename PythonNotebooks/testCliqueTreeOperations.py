import sys
import numpy as np
from Factor import *
from PGMcommon import *
from PedigreeFactors import *
from FactorOperations import *
from GeneticNetworkFactory import *
from CliqueTree import *
import networkx as nx
import matplotlib.pyplot as plt
from CliqueTreeOperations import *
#to create a clique tree, we start with a list of factors
#and potentially some observed evidence
alphaList=[.8,.6,.1]
allelefreq=[.1,.9]
chrom='12'
position=1000
g1=GeneticNetworkFactory('sixperson.ped',alphaList,allelefreq, chrom,position)
g1.constructNetwork()
factorList=g1.getFactorList()

cTree = createCliqueTree(factorList)

G=nx.from_numpy_matrix( cTree.getEdges() )
nx.draw_shell(G)
plt.show()

prunedCTree=PruneTree( cTree )
G=nx.from_numpy_matrix( prunedCTree.getEdges() )
nx.draw_shell(G)
plt.show()
