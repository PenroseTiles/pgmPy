from Factor import *
from PGMcommon import *
from CliqueTree import *
from CliqueTreeOperations import *
from FactorOperations import *
import scipy.io as sio
import numpy as np
import pprint
import pdb

matfile='/Users/indapa/software/PGM/Prog4/PA4Sample.mat'
mat_contents=sio.loadmat(matfile)

mat_struct=mat_contents['SixPersonPedigree']
val=mat_struct[0]

factorList=[]
for elem in val:
    
    (var, card, val) =elem
    f= Factor( var[0].tolist(), card[0].tolist(), val[0].tolist(), 'factor' )
    #print f
    factorList.append(f)
#print

jointFactor = ComputeJointDistribution(factorList)
#print jointFactor

MARGINALS= ComputeExactMarginalsBP( factorList, [], 1)
P=CreatePrunedInitCtree(factorList, [] )
(P, MESSAGES) = CliqueTreeCalibrate(P, isMax=1)
jointDistribution=ComputeJointDistributionFromCalibratedCliqueTree(P, MESSAGES, 1)

print jointDistribution


    #for marginal in MARGINALS:
    #print marginal
#print

MAPAssignment=MaxDecoding( MARGINALS  )
print MAPAssignment


    #for m in MARGINALS:
    #m.setVal( np.log( lognormalize(m.getVal()   )   ) )
#print np.sum( lognormalize(m.getVal() ) )

#print MAPAssignment

    #for marginal in MARGINALS:
    #print marginal.getVal()
#print


