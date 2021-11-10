import time as time
import numpy as np
import random as random


"""
The module containing functions to implement MCMC methods for Ising models.
"""
def update_external(current,beta,A,b,current_energy):
    """
    Implements one step of the Metropolis Markov chain algorithm to sample from a Gibbs state of a classical Ising model.

    Args:
        current (array with +/-1 entries): initial state of the chain.
        beta(float): inverse temperature of the Gibbs state we wish to sample from.
        A(scipy sparse matrix): the (weighted) adjacency matrix of the Ising model whose energy we wish to minimize.
        b(vector): the external fields.
        currrent_energy (float): energy of the initial state current with respect to the model.

    Returns:
        [current,new_energy] ([array with +/-1 entries, float]): a new candidate state for the spin system and its corresponding energy.


        
    """
    
    n=A.shape[0]
    new_energy=current_energy
    x=random.randint(0,n-1)  
    sum_neigh=0
    B=A[x,:].nonzero()
    for k in range(0,len(B[1])):
        sum_neigh+=A[x,B[1][k]]*current[B[1][k]]
    
    change_energy=-4*current[x]*sum_neigh-4*current[x]*b[x]
    #print(current,x,current[x],2*current[x]*sum_neigh)
    if change_energy<0:
        current[x]=(-1)*current[x]
        new_energy=current_energy+change_energy
    else:
        p=random.random()
        if np.exp(-beta*(change_energy))>p:
            current[x]=(-1)*current[x]
            new_energy=current_energy+change_energy
    return [current,new_energy]


def telescopic_product_external(A,b,schedule,samples,burn_in,verbose=1):
    """
    Uses the telescopic product approach to compute the partition function of an Ising model for a list of inverse temperatures. See the page of the Monte Carlo methods for an explanation of what this means exactly.



    Args:
        A(scipy sparse matrix): the (weighted) adjacency matrix of the Ising model whose energy we wish to minimize.
        b(vector): the external fields.
        schedule(array of floats): an array containing the annealing schedule. That is, a list of inverse temperatures whose partition function we wish to estimate.
        samples(integer): the number of samples that we will take for the estimation at every value of the inverse temperature.
        burn_in(integer): how many steps of the Metropolis Markov chain we will take to converge approximately to the Gibbs state distribution and start estimating the ratio of the partition function. 

    
    
    Returns:
        estimates_ratio ([array of floats]): array containing the ratios of the partition functions in the schedule. That is, we compute :math:`\\frac{Z_{\\beta_{i+1}}}{Z_{\\beta_{i}}}` for consecutive entries :math:`\\beta_i,\\beta_{i+1}` of the annealing schedule.

        
    """
    estimates_ratio=[]
    n=A.shape[0]
    first=1
    for ratio in range(0,len(schedule)-1):
        
        if ratio==1 and verbose==1:
            print("It will take roughly ",int((end-start)*len(schedule)/60),"minutes to complete the estimates")
        start = time.time()
        current_ratio=0
        for m in range(0,samples):
            
            
            
            if first==1:
                current=2*np.random.randint(2,size=n)-np.ones(n)
                blob=A.dot(current)
                energy=current.dot(blob)+current.dot(b)
            if first<1:
                burn_in=5
            first=0
            for  k  in range(0,burn_in):
                [current,energy]=update_external(current,schedule[ratio],A,b,energy)
            end=time.time()
            
            current_ratio+=np.exp(-(schedule[ratio+1]-schedule[ratio])*energy)
        estimates_ratio.append(current_ratio/samples)
        end=time.time()
        #print("it took:",end-start)
            
    return estimates_ratio
