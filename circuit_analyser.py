
import numpy as np
import warnings

from numpy import linalg as LA
from scipy import linalg
import cirq
"""
The module containing functions to estimate the relative entropy.
""" 

def entropy_output(circuit,number_qubits,fixed_point,contraction):
    """
    Estimates the relative entropy between the output of a noisy circuit and the fixed point of the noise. The code assumes that the initial state of the circuit is :math:`\\ket{0}^{\\otimes n}`. Note that all arguments related to the noise are automatically provided by the quantum channel class.

    Returns:
        final_ent (float): upper bound on the relative entropy between output of the circuit and fixed point of the noise.

    Args:
    
        circuit(Cirq circuit): cirq noiseless circuit whose output's energy we wish to lower-bound.
        number_qubits (integer): the number of qubits of the circuit.
        fixed_point (2x2 complex array): density matrix corresponding to the fixed point of the noise.
        contraction (float): contraction of the noise towards the fixed point. 
    """
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
                
            
    final_ent=initial*(contraction**(number_moments))+total_contrib_moments
    return final_ent







