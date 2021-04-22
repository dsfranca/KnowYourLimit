#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 22:34:45 2021

@author: danielstilckfranca
"""
import numpy as np
import warnings

from numpy import linalg as LA
from scipy import linalg
import cirq
#computes the relative entropy between fixed point of noise and output given the fixed point, contraction rate
#and cirq circuit. Will assume that the circuit is initiated at 0. 
def entropy_output(circuit,number_qubits,fixed_point,contraction):
    #initial relative entropy, given by n*tr(|0><0|log(fixed_point))
    initial=-number_qubits*linalg.logm(fixed_point)[0,0]
    #print(initial)
    #this part will record by how much the circuit leaves the fixed point invariant at each stage
    unitary_part=[]
    number_moments=0
    #square root of inverse to compute max entropy
    square_root_inv=linalg.sqrtm(linalg.inv(fixed_point))
    #2 copies of square root
    square_root_inv2=np.kron(square_root_inv,square_root_inv)
    #2 copies of usual
    fixed_point2=np.kron(fixed_point,fixed_point)
    
    contrib_moments=[]
    for moment in circuit:
        number_moments+=1
        contrib_current_moment=0
        for gate in moment:
            
            U=np.matrix(cirq.unitary(gate))
            if np.shape(U)[0]==4:
                diff=square_root_inv2@U@fixed_point2@U.getH()@square_root_inv2
                contrib_current_moment+=np.log(linalg.norm(diff,2))
            else:
                diff=square_root_inv@U@fixed_point@U.getH()@square_root_inv
                contrib_current_moment+=np.log(linalg.norm(diff,2))
        contrib_moments.append(contrib_current_moment)
        
    total_contrib_moments=0
    for k in range(0,number_moments):
        #print("new contribution",contraction**(number_moments-k)*contrib_moments[k])
        total_contrib_moments+=contraction**(number_moments-k)*contrib_moments[k]
                
            
        
    return initial*(contraction**(number_moments))+total_contrib_moments







