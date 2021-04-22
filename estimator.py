#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 10:43:35 2021

@author: danielstilckfranca
"""
from circuit_compiler import *
from partition_function_estimator import *
from noise_model import *
from circuit_analyser import *
import cirq
import warnings
import cvxgraphalgs as cvxgr
from cvxgraphalgs.structures import Cut
from scipy import sparse

#lower bounds the expectation value of the output of a compiled and routed QAOA circuit 
#by computing the  partition funciton brute force.
def estimator_brute_force(noise_model,problem_graph,device,gammas,betas):
    
    noise=quantum_channel(noise_model)
    sigma=noise.fixed_point
    circuit=compiled_routed_qaoa(problem_graph,gammas,betas,device)
    n_qubits=problem_graph.order()
    if n_qubits>15:
        warnings.warn("You are trying to compute the partition function of more than 15 qubits. It may take substantial time")
    
    
    rel_ent=entropy_output(circuit,n_qubits,sigma,noise.contraction)
    A=nx.adjacency_matrix(problem_graph)
    output_energy=np.real(estimate_output_energy_brute(A,sigma,rel_ent,6))
    ground_state_energy=ground_state_brute_force(A)
    print(output_energy,ground_state_energy,noise.contraction)
    
    
    
    return output_energy,ground_state_energy

def estimator_brute_force_weighted(noise_model,problem_graph,device,gammas,betas):
    
    noise=quantum_channel(noise_model)
    sigma=noise.fixed_point
    circuit=compiled_routed_weighted_qaoa(problem_graph,gammas,betas,device)
    n_qubits=problem_graph.order()
    if n_qubits>15:
        warnings.warn("You are trying to compute the partition function of more than 15 qubits. It may take substantial time")
    
    
    rel_ent=entropy_output(circuit,n_qubits,sigma,noise.contraction)
    A=nx.adjacency_matrix(problem_graph)
    #change later, done to account for measurement error
    output_energy=np.real(estimate_output_energy_brute(A,sigma,0.9*rel_ent,6))
    ground_state_energy=ground_state_brute_force(A)
    
    
    
    return output_energy,ground_state_energy



def estimator_monte_carlo(noise_model,problem_graph,device,gammas,betas,verbose=1):
    
    noise=quantum_channel(noise_model)
    sigma=noise.fixed_point
    if verbose==1:
        print("Compiling circuit")
    circuit=compiled_routed_weighted_qaoa(problem_graph,gammas,betas,device)
    if verbose==1:
        print("Circuit compiled")
    n_qubits=problem_graph.order()
    A=nx.adjacency_matrix(problem_graph)
    A=A.asfptype()
    [v,s,w]=sparse.linalg.svds(A,k=1)
    normA=s[0]
    
    #determine external field
    sigma=noise.fixed_point
    y=sigma[0,0]-sigma[1,1]
    #weird normalization. I have encountered some versions of numpy that use base 2.
    beta_sigma=(0.5/np.log(np.exp(1)))*np.log((1+y)/(1-y))
    b=0.25*beta_sigma*np.ones(n_qubits)
    if verbose==1:
        print("Estimating relative entropy of output")
    rel_ent=entropy_output(circuit,n_qubits,sigma,noise.contraction)
    if verbose==1:
        print("Estimating partition functions")
    schedule=np.linspace(rel_ent/(n_qubits*10),5/normA,2*n_qubits)
    partition=telescopic_product_external(A,b,schedule,50*n_qubits,8*n_qubits)
    output_energy=np.real(estimate_output_energy_monte_carlo(schedule,partition,rel_ent))
    #ground_state_energy=ground_state_brute_force(A)
    #print(output_energy,ground_state_energy,noise.contraction)
    
    
    
    return output_energy #output_energy,ground_state_energy





def SDP_average(problem_graph,samples):
    results=[]
    for k in range(0,samples):
        sdp_cut = cvxgr.algorithms.goemans_williamson_weighted(problem_graph)
        A=nx.adjacency_matrix(problem_graph)
        results.append(2*(0.5*np.sum(A)-2*sdp_cut.evaluate_cut_size(problem_graph)))
    return np.sum(results)/samples
    
#    
#noise_model = cirq.asymmetric_depolarize(
#    p_x=0.01,
#    p_y=0.01,
#    p_z=0.01,
#)
#noise_model=cirq.channel(noise_model)
#m=4
#gammas=[0.1,0.5,0.7,0.9]
#betas=[0.7,0.5,0.1,0.05]
#tk_device=pytket.device.Device({}, {}, pytket.routing.SquareGrid(10,10))  
#problem_graph = nx.random_regular_graph(d=4, n=100)
#
#
##[estimate,true_value]=estimator_brute_force(noise_model,problem_graph,tk_device,gammas,betas)
#hope,betas,partition=estimator_monte_carlo(noise_model,problem_graph,tk_device,gammas,betas)
#print("Estimate",hope)
#
#A=nx.adjacency_matrix(problem_graph)
#print("Exact",partition_brute_force(A,np.eye(2)/2,schedule[-1]))
#
###define noise_model  which is mixture of dep and amplitude damping
#noise_model2=cirq.amplitude_damp(0.01)
#
#noise_model2=cirq.channel(noise_model2)
#new_noise=[]
#for A in noise_model2:
#    new_noise.append(np.sqrt(0.5)*A)
#for A in noise_model:
#    new_noise.append(np.sqrt(0.5)*A)
#
##[estimate,true_value]=estimator_brute_force(new_noise,problem_graph,tk_device,gammas,betas)  
#
#noise_model=quantum_channel(new_noise)
#
#[schedule,partition]=estimator_monte_carlo(new_noise,problem_graph,tk_device,gammas,betas)
#print("Estimate",-np.sum(np.log(partition)))
#
#A=nx.adjacency_matrix(problem_graph)
#
#print("Exact",partition_brute_force(A,np.real(noise_model.fixed_point),schedule[-1]))


#
#
#
##
##
#sdp_cut = cvxgr.algorithms.goemans_williamson_weighted(problem_graph)
#
#
#nx.adjacency_matrix(problem_graph)
#
#
#x=SDP_average(problem_graph,1)
##print("Exepcted sdp",x)
#print("Approximation ration",x/true_value)



#
#
#
#m=4
#n_system=m**2
#problem_graph = nx.random_regular_graph(d=4, n=n_system)
#    
#tk_device=pytket.device.Device({}, {}, pytket.routing.SquareGrid(m,m))  
#gammas=np.ones(2)
#betas=np.ones(2)
#circuit=compiled_routed_qaoa(problem_graph,gammas,betas,tk_device)
#    
#
#problem_graph=nx.complete_graph(n_system)
#
#
#gammas=np.ones(2)
#betas=np.ones(2)
#circuit=compiled_routed_qaoa(problem_graph,gammas,betas,tk_device)
#
#
#
#
#n_system=12
#problem_graph = nx.random_regular_graph(d=2, n=n_system)
#A=nx.adjacency_matrix(problem_graph)
#sigma=np.array([[0.5,0],[0,0.5]])
#
#
#print(energy_sigma(A,sigma,10000))
#
#beta=1
#print("Ground state energy",partition_brute_force(A,sigma,3)/3)
#print("Max mixed",estimate_output_energy_brute(A,sigma,n_system*0.1,3))
#
#sigma=np.array([[0.99,0],[0,0.01]])
#
#
#print(energy_sigma(A,sigma,10000))
#limit=-np.ones(n_system)
#print(limit@A@limit)
#
#beta=0.1
#print("Amp damping",estimate_output_energy_brute(A,sigma,n_system*0.01,3))