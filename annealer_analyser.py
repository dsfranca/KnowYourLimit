import scipy.integrate as integrate
import scipy.special as special
import numpy as np

def entropy_output_annealer(schedule,T,number_qubits,fixed_point,contraction):
    """
    Estimates the relative entropy between the output of a noisy annealer and the fixed point of the noise. The code assumes that the initial state of the circuit is :math:`|+\\rangle^{\\otimes n}` and that the annealing schedule is of the form :math:`f(t/T)H_0+g(t/T)H_I`, where :math:`f` is the annealing schedule and :math:`H_0=-\\sum_iX_i`. The function :math:`g` does not contribute to the final entropy, so we ignore if for now. Note that the contraction is now the contraction rate, as we are in continuous time.

    Returns:
        final_ent (float): upper bound on the relative entropy between output of the circuit and fixed point of the noise.

    Args:
    
        schedule(callable function): function :math:`f:[0,1]\\to\mathbb{R}` that describes the annealing schedule.
        T(float): total annealing time.
        number_qubits (integer): the number of qubits of the annealer.
        fixed_point (2x2 complex array): density matrix corresponding to the fixed point of the noise. Note that is has to correspond to a classical state and we adopt the convention that the probability of observing 0 is larger than 1.
        contraction (float): contraction rate of the Linfbladian towards the fixed point. Note that we have 
    """
    #determine inverse temperature of fixed point
    y=fixed_point[0,0]
    gamma=0.5*np.log(y/(1-y))
    #initial relative enropy between fixed point and |+>.
    initial_rel_ent=number_qubits*np.log(2*np.cosh(gamma))
    #commutator of the fixed point with the sigma_x part of the Hamiltonian.
    commutator=2*np.sinh(gamma)
    #compute the relative entropy decay term
    rel_ent_decay=initial_rel_ent*np.exp(-contraction*T)
    #compute the relative entropy decay from the commutator
    result_commutator = commutator*number_qubits*integrate.quad(lambda x: schedule(x/T)*np.exp(-contraction*(x-T)), 0,T)[0]
    final_ent=rel_ent_decay+result_commutator
    return final_ent


    

