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
mat_struct=mat_contents['FactorMax']
val=mat_struct[0,0]
input_factors = val['INPUT1'][0][0]
var = input_factors[0].flatten().tolist()
card=input_factors[1].flatten().tolist()
value=input_factors[2].flatten().tolist()
print var
print card
print value
INPUT1= Factor( var, card, value, 'test')
INPUT2= val['INPUT2'].flatten()
print INPUT1
print INPUT2
print FactorMaxMarginalization(INPUT1, INPUT2)
#example used in section 13.2 pg 555 of Friedman and Koller
print "====="
psi=Factor( [ 1,2,3], [3,2,2], [.25,.05,.15,.08,0,.09,.35,.07,.21,.16,0,.18])
maxfactor= FactorMaxMarginalization(psi, [2])
print maxfactor
print IndexToAssignment(np.arange(6),[3,2])
