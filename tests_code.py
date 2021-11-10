
import warnings
import numpy as np
import quimb.tensor as qtn
#from quimb.tensor.tensor_arbgeom import TensorNetworkGen, TensorNetworkGenVector

from quimb.utils import concat
#import random
import scipy
#import time
#import itertools
from scipy.sparse import csr_matrix, find
#import random
import numpy as np
#import time
#import itertools
#import pandas as pd

#import math
#import collections
#import cirq

#import cvxgraphalgs as cvxgr
#from cvxgraphalgs.structures import Cut


from partfunc_estimators import *
#compute he partiion funcion of product state.
n=10
import networkx as nx
G=nx.random_regular_graph(d=4,n=n)

A=0.0001*nx.to_scipy_sparse_matrix(G)
bs=np.ones(n)
beta0=0
beta1=1
[betas,parts]=partition_function_estimator_Ising(beta0,beta1,A,'TN',bs,step_size=0.1)

print(betas)
print(parts)
print(n*np.log(np.exp(betas)+np.exp(-betas)))