import scipy.io as sio
import numpy as np
import pprint
from Factor import *
from PGMcommon import *
from CliqueTree import *
from CliqueTreeOperations import *
from FactorOperations import *
import pickle

""" load the test data from the matlab file """

matfile='/Users/amit/BC_Classes/PGM/Prog4/PA4Sample.mat'
mat_contents=sio.loadmat(matfile)
mat_struct=mat_contents['GetNextC']
np.shape(mat_struct)
val=mat_struct[0,0] 

input_edges = val['INPUT1']['edges'][0][0]
input_cliqueList= val['INPUT1']['cliqueList'][0][0][0]
clique_list_factorObj=[]
for tpl in input_cliqueList:
    (var, card, values)=tpl
    f= Factor( var[0].tolist(), card[0].tolist(), values[0].tolist(), 'factor' )
    clique_list_factorObj.append(f)

messages=val['INPUT2']
messageFactors=[]
(nrow,ncol)=np.shape(messages)
for i in range( nrow ):
    for j in range ( ncol ):
        (var, card, values)=messages[i][j]
        f=Factor( var[0].tolist(), card[0].tolist(), values[0].tolist(), 'factor' )
        messageFactors.append(f)
MESSAGES= np.reshape( np.array ( messageFactors ), (nrow,ncol) )


""" dump it to Python pickle file """
with open('GetNextC.INPUT1.pickle', 'wb') as f:
    pickle.dump(clique_list_factorObj,f)
with open('GetNextC.INPUT2.pickle', 'wb') as f:
    pickle.dump(MESSAGES, f)

P=CliqueTree( clique_list_factorObj ,  input_edges, clique_list_factorObj, [])
(a,b)=getNextClique(P,MESSAGES)
print 'a: ', a, ' b:', b





