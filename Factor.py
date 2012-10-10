import sys
import numpy as np
class Factor(object):
    

    def __init__(self, var=[], card=[], val=[], name= 'None'):
        self.var= np.array(var)
        self.card=np.array(card)
        self.val=np.array(val)
        self.name=name

    def __str__(self):
        varstring= " ".join ( map(str, self.var) )
        cardstring=" ".join ( map(str, self.card) )
        valstring= " ".join( map(str, self.val))
        return "\n".join( [ self.name,varstring, cardstring, valstring])


    def setVar(self, var):
        self.var=np.array(var)

    def getVar(self):
        return self.var

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




