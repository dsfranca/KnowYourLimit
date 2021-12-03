import scipy.integrate as integrate
import scipy.special as special
import numpy as np

def entropy_output_annealer(schedule,T,number_qubits,fixed_point,contraction):
    """
    Estimates the relative entropy between the output of a noisy annealer and the fixed point of the noise. The code assumes that the initial state of the circuit is :math:`\\ket{+}^{\\otimes n}` and that the annealing schedule is of the form :math:`f(t/T)H_0+g(t/T)H_I`, where f is the annealing schedule and :math:`f(t/T)H_0=-\sum_iX_i'. The function :math:'g' does not contribute to the final entropy, so we ignore if for now. Note that the contraction is now the contraction rate, as we are in continuous time.

    Returns:
        final_ent (float): upper bound on the relative entropy between output of the circuit and fixed point of the noise.

    Args:
    
        schedule(callable function): function :math:`f:[0,1]\\to\mathbb{R}` that describes the annealing schedule.
        T(float): total annealing time.
        number_qubits (integer): the number of qubits of the annealer.
        fixed_point (2x2 complex array): density matrix corresponding to the fixed point of the noise.
        contraction (float): contraction rate of the Linfbladian towards the fixed point. Note that we have 


