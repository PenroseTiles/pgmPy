from Factor import *
from PGMcommmon import *
from PedigreeFactors import *
import numpy as np


class CliqueTree(object):
    """ represent a CliqueTree (Junction Tree)
        a clique tree is a graph whose nodes are cliques
        has an adjacency matrix and a factorList from which
        the graph was built from
        a clique is associated with a subset of variables
        an edge connects to cliques is associated with the
        sepset (intersection of variables between) the  cliques"""

    def __init__(self, nodeList=[], edges=[], factorList=[]):
        self.nodeList=nodeList
        self.edges=edges
        self.factorList=factorList

    def setNodeList(self,nodeList):
        self.nodeList=nodeList
    def setEdges(self,edges):
        self.edges=edges
    def setFactorList(self,factorList):
        self.factorList=factorList

    def getNodeList(self):
        return self.nodeList
    def getEdges(self):
        return self.edges
    def getFactorList(self):
        return self.factorList



