#!/usr/bin/env python
from Factor import *
import numpy as np

def IndexToAssignment( I, D):
    a=np.reshape ( np.arange(np.prod(D)).repeat(len(D)), (np.prod(D),len(D)))
    #print np.shape(a)
    print a

    b=tmp=list( D[:-1] )
    tmp.insert(0,1)
    tmp =np.cumprod ( np.array (tmp) )
    print tmp
