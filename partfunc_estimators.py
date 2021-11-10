import networkx as nx
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

"""
The module containing functions that estimate partition functions of Ising models.
"""
class weights_graph:
    def __init__(self,G):
        self.G=G
    def __call__(self,k,l):
        return -self.G[k][l]['weight']
class external_field:
    def __init__(self,b):
        self.b=b
    def __call__(self,k):
        return self.b[k]

def partition_function_estimator_Ising(beta0,beta1,A,method,b=0,step_size=0.05):
    """
    Returns the value of the log partition function of the Ising model defined by the matrix A and with external field b for all inverse
    temperatures in the interval [beta0,beta1]. The partition function is evaluated at every step-size.
    The method can either be Tensor networks or Monte Carlo. For Monte Carlo methods, the annealing schedule and the number of samples will be set automatically.
    Note that for Monte Carlo methods it could be the case that the estimate is not reliable for high values of beta1. Tensor network methods, on the other hand,
    cannot handle very large instances.

    Returns:
        [betas,partitions] ([list of floats, list of floats]): a list of inverse temperatures and the corresponding log partition functions.

    Args:
    
        beta0(float): initial inverse temperature.
        beta1(float): final inverse temperature.
        A(scipy sparse matrix): the (weighted) adjacency matrix of the Ising model whose energy we wish to minimize.
        b(vector): the external fields.
        method(str): either 'TN' or 'MC'. Determines which method will be used to compute the partition function. TN is tensor network, MC is Monte Carlo.
        step_size(float): interval at which we evaluate the partition function.
        
    """
    n=A.shape[0]
    G=nx.from_scipy_sparse_matrix(A)
    if method is 'TN' and n>75:
        warnings.warn("You are conctracting a tensor network with more than 75 nodes. This might take long")
    
    if method is 'TN':
        steps=int((beta1-beta0)/step_size)
        betas=np.linspace(beta0,beta1,steps)
        partitions=[]
        weight_func=weights_graph(G)
        if b is not 0:
            external=external_field(b)
        for beta in betas:
            if b is not 0:
                tn=qtn.TN_classical_partition_function_from_edges(G.edges(), beta=beta,j=weight_func,h=external)
                partitions.append(np.log(tn.contract()))
    
            else: 
                tn=qtn.TN_classical_partition_function_from_edges(G.edges(), beta=beta,j=weight_func)
                partitions.append(np.log(tn.contract()))
    return [betas,np.real(partitions)]

