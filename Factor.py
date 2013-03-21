import sys
import numpy as np
class Factor(object):

    """ Represents a factor in a PGM. A factor has variable scope, cardinality, value, and potentially a name.
        A factor's values (which can potentially be multi-dimensional table, are represented as NumPy 1d-arrays.
        See Chapter10 Friedman pg. 358 and Koller for more info. """



    def __init__(self, var=[], card=[], val=[], name= 'None'):
        """ a factor has list of variables, each with a cardinality, and for each possible assignment to its variable(s),
        a position in the val array."""
        self.var= np.array(var)
        self.card=np.array(card)
        self.val=np.array(val)
        self.name=name

    def __str__(self):
        varstring= " ".join ( map(str, self.var) )
        cardstring=" ".join ( map(str, self.card) )
        valstring= " ".join( map(str, self.val))
        return "\n".join( [ 'name: ' + self.name,'var: '+ varstring, 'card: '+ cardstring, 'val: ' + valstring])


    def setVar(self, var):
        self.var=np.array(var)

    def getVar(self):
        return self.var

    def getVarCount(self):
        return len( self.var.tolist() )

    def setVal(self, val):
        self.val=np.array(val)

    def getVal(self):
        return self.val

    def setCard(self,card):
        self.card=np.array(card)

    def getCard(self):
        return self.card

    def getName(self):
        return self.name

    def setName(self, name):
        self.name=name




