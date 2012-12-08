from Factor import *
from PGMcommon import *
from CliqueTree import *
from CliqueTreeOperations import *
from FactorOperations import *
import scipy.io as sio
import numpy as np
import pprint
import pdb
matfile='/Users/amit/BC_Classes/PGM/Prog4/PA4Sample.mat'
mat_contents=sio.loadmat(matfile)
mat_struct=mat_contents['ExactMarginal']
val=mat_struct[0,0]
input_factors = val['INPUT'][0]
factorList=[]
for tpl in input_factors:
    (var, card, values)=tpl
    #print var, card, values
    f= Factor( var[0].tolist(), card[0].tolist(), values[0].tolist(), 'factor' )
    #print f
    factorList.append( f )

ComputeExactMarginalsBP( factorList )
