from Factor import *
from PGMcommon import *
from CliqueTree import *
from CliqueTreeOperations import *
from FactorOperations import *
import scipy.io as sio
import numpy as np
import pprint
import pdb
import matplotlib.pyplot as plt
import networkx as nx
matfile='/Users/amit/BC_Classes/PGM/Prog4/PA4Sample.mat'
mat_contents=sio.loadmat(matfile)
mat_struct=mat_contents['OCRNetworkToRun']
val=mat_struct
factorList=[]

ALPHABET=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

for data in val:
    
    (var, card, values)=data[0]
    f= Factor( var[0].tolist(), card[0].tolist(), values[0].tolist(), 'factor' )
    factorList.append( f )

#MARGINALS= ComputeExactMarginalsBP( factorList, [], 1 )
#MAPAssignment=MaxDecoding( MARGINALS )
#print "".join( [ALPHABET[idx] for idx in MAPAssignment] )

MARGINALS= ComputeExactMarginalsBP( factorList, [], 1 )
for m in MARGINALS:
    log_val= m.getVal()
    prob_val_normalized=np.log( lognormalize( log_val ) )
    m.setVal(prob_val_normalized)


MAPAssignment=MaxDecoding( MARGINALS  )
print "".join( [ALPHABET[idx] for idx in MAPAssignment] )

for m in MARGINALS:
    print np.sum( lognormalize(m.getVal() ) )


#V=getUniqueVar(factorList)
#print 'unique variables:'
#print V

#cTree=CreatePrunedInitCtree(factorList)
#G=nx.from_numpy_matrix( cTree.getEdges() )
#nx.draw_shell(G)
#plt.show()
