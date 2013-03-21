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
mat_struct=mat_contents['MaxDecoded']
val=mat_struct[0,0]
input_factors = val['INPUT']
factorList=[]



for elm in input_factors:
    #print ary[0]
    #print
    (var, card, values)=elm[0]
    #print var, card, values
    f= Factor( var[0].tolist(), card[0].tolist(), values[0].tolist(), 'factor' )
    print f
    factorList.append( f )

DECODE= MaxDecoding( factorList )
ALPHABET=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

print DECODE
print [ALPHABET[idx] for idx in DECODE]
