class CliqueTree(object):
    'MyClass constructor'

    def __init__(self, nodeList=[], edges=[], factorList=[],evidence=[]):
        self.nodeList=nodeList
        self.edges=edges
        self.factorList=factorList
        self.evidence=evidence

    def setNodeList(self,nodeList):
        self.nodeList=nodeList

    def setEdges(self,edges):
        self.edges=edges

    def setFactorList(self,factorList):
        self.factorList=factorList

    def setEvidence(self,evidence):
        self.evidence=evidence

    def getNodeList(self):
        return self.nodeList

    def getEdges(self):
        return self.edges

    def getFactorList(self):
        return self.factorList

    def getEvidence(self):
        return self.evidence

    def eliminateVar(self,Z):
        """ a variable elimination function
            based on https://github.com/indapa/PGM/blob/master/Prog4/EliminateVar.m """

        useFactors = []#the index of the factor that contains the variable Z
        scope = []

        for i in range (len(self.factorList)):
            if Z in self.factorList[i].getVar().tolist():
                useFactors.append(i)
                scope=list(set.union(set(scope), self.factorList[i].getVar().tolist() ))

        print useFactors
        print scope

        # update edge map
        # These represent the induced edges for the VE graph.

        for i in range ( len(scope)):
            for j in range ( len(scope)):
                if i!=j:
                    self.edges[i,j]=1
                    self.edges[j,i]=1
        self.edges[Z,:]=0
        self.edges[:,Z]=0
        

