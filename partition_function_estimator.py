#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 09:00:01 2021

@author: danielstilckfranca
"""

from scipy.sparse import csr_matrix, find
import random
import numpy as np
import time
import itertools
import pandas as pd
import networkx as nx

def partition_brute_force(A,sigma,beta):
    n=A.shape[0]
    estimate=0
    #determine inverse temperature of fixed point, assuming it is classical
    y=sigma[0,0]-sigma[1,1]
    #weird normalization. I have encountered some versions of numpy that use base 2.
    beta_sigma=(0.5/np.log(np.exp(1)))*np.log((1+y)/(1-y))
    #print("Invers temperature",beta_sigma)
    for initial in itertools.product([-1, 1], repeat=n):
        initial=np.array(initial)
        external_field_contrib=np.sum(initial)
        #add part stemming from A
        new_energy=A.dot(initial)
        #add part stemming from external field
        
        #print(new_energy.shape,initial.shape)
        new_energy=initial.dot(new_energy)
        
        estimate=estimate+np.exp(-beta*new_energy+beta_sigma*external_field_contrib)
        #print(initial,new_energy,estimate/(2**n))
     
    #print("I believe this:",estimate)
        
    return n*np.log((2*np.cosh(beta_sigma)))-np.log(estimate)

#computes energy at limit
def energy_sigma(A,sigma,samples):
    estimates=[]
    n=A.shape[0]
    state=np.ones(n)
    for sample in range(0,samples):
        for k in range(0,n):
            p=random.random()
            if p<sigma[0,0]:
                state[k]=-1
            else:
                state[k]=1
        
        new_energy=A.dot(state)
        new_nergy=state.dot(new_energy)
        estimates.append(new_nergy)
    return np.sum(estimates)/samples

def estimate_output_energy_brute(A,sigma,rel_ent,beta_max):
    betas=np.linspace(rel_ent/10,beta_max,100)
    estimate=[]
    for beta in betas:
        part=partition_brute_force(A,sigma,beta)
        estimate.append((part-rel_ent)/beta)
    return max(estimate)

def estimate_output_energy_brute(A,sigma,rel_ent,beta_max):
    betas=np.linspace(rel_ent/10,beta_max,100)
    estimate=[]
    for beta in betas:
        part=partition_brute_force(A,sigma,beta)
        estimate.append((part-rel_ent)/beta)
    return max(estimate)


def estimate_output_energy_monte_carlo(betas,partitions,rel_ent):
    
    estimate=[]
    for k in range(0,len(betas)):
        part=-np.sum(np.log(partitions[0:k]))
        estimate.append((part-rel_ent)/betas[k])
    return max(estimate)    

def ground_state_brute_force(A):
    n=A.shape[0]
    estimate=0

    
    for initial in itertools.product([-1, 1], repeat=n):
        initial=np.array(initial)
        external_field_contrib=np.sum(initial)
        #add part stemming from A
        new_energy=A.dot(initial)
        #add part stemming from external field
        
        #print(new_energy.shape,initial.shape)
        new_energy=initial.dot(new_energy)
        
        if new_energy<estimate:
            estimate=new_energy
        #print(initial,new_energy,estimate/(2**n))
     
        
    return estimate

def update_external(current,beta,A,b,current_energy):
    
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
            
            
def telescopic_product_external(A,b,schedule,samples,burn_in):
    estimates_ratio=[]
    n=A.shape[0]
    first=1
    for ratio in range(0,len(schedule)-1):
        
        if ratio==1:
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









