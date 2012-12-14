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
mat_struct=mat_contents['MaxSumCalibrate']
val=mat_struct[0,0]
input_edges = val['INPUT']['edges'][0][0]
#print input_edges
input_cliqueList= val['INPUT']['cliqueList'][0][0][0]
clique_list_factorObj=[]
for tpl in input_cliqueList:
    (var, card, values)=tpl
    f= Factor( var[0].tolist(), card[0].tolist(), values[0].tolist(), 'factor' )
    #print f 
    clique_list_factorObj.append(f)

P=CliqueTree( clique_list_factorObj ,  input_edges, clique_list_factorObj, [])

P=CliqueTreeCalibrate(P,1)

for f in P.getNodeList():
    print f
    print "=="
