################
Know Your Limits
################


What is KnowYourLimits?
------------------------
KnowYourLimits is a Python library to benchmark the performance of various quantum algorithms under noise. It is based on the methods developed in the paper
`"Limitations of optimization algorithms on noisy quantum devices" <https://doi.org/10.1038/s41567-021-01356-3>`_.
The intent of this package is to make the tools developed in that work easily accessible to the quantum computation community and people interested in estimating wether 
near term quantum devices have the potential to provide better solutions for their applications.

This first version of the project will focus on quantum algorithms that have as a goal to minimize the energy of an `Ising model <https://en.wikipedia.org/wiki/Ising_model>`_.
This is justtified by several reasons:

* The polyvalence of the Ising model: it is well-known that many different problems in combinatorial optimization and physics can be formulated as minimizing the energy of an Ising model. See e.g. `this paper <https://www.frontiersin.org/articles/10.3389/fphy.2014.00005/full>`_ for examples.
* Many proposals for the use of near-term devices have this as their goal to solve this problem. Some prominent examples that our package can handle are the Quantum Approximate Optimization Algorithm (see e.g. `this page <https://marwahaha.github.io/qaoa-reference/>`_ for references) and adiabatic quantum computing (see e.g. `this review <https://arxiv.org/pdf/1611.04471.pdf>`_).
* This is the setting in which the techniques of "Limitations of optimization algorithms on noisy quantum devices" are most effective.
  
The performance of near-term quantum devices on solving optimization problems is limited by the interplay of two factors:

#. Their inherent noise and lack of error correction.
#. Their limited connectivity. That is, it is not possible to apply gates between all possible qubits of the system.

Point 1. cleary limits how many operations we can perform before the output of the noisy quantum device is dominate by the noise, whereas point 2. amplifies this effect.
Indeed, as the limited connectivity then requires us to apply extra operations to ensure that we implement a given quantum circuit. Thus, KnowYourLimits needs to 
take points 1. and 2 into acccount.


What does KnowYourLimits do?
----------------------------
KnowYourLimits mostly uses `Cirq package <https://quantumai.google/cirq>`_ objects to specify noise models and quantum circuits.

KnowYourLimits takes the following inputs:

#. A noise model for the underlying device.
#. A graph of the connectivity of the device.
#. The circuit to be implemented or an annealing schedule.
#. The Hamiltonian of the Ising state whose energy we wish to minimize.

Given this input, it uses `pytket <https://github.com/CQCL/pytket>`_ to route the circuit the device's architecture. From that it then outputs a **lower bound to the 
energy of the output** of the noisy circuit. KnowYourLimits can then also compare this lower bound with the output of simple classical algorithms for the problem at hand.
This way the user can infer whether there is room for better solutions by using that given noisy quantum device or classical methods will still outperform noisy quantum devices.

For some noise models the exact circuit being implemented does not influence the lower bound, only its depth. This is the case for noise models that drive the 
device to the maximally mixed state, like depolarizing noise.

Importantly, the method used by KnowYourLimits does not require simulating the underlying quantum circuit. This way it can still give estimates for quantum devices 
comprised of up to thousands of qubits within a reasonable time, while simulating the quantum devices is out of reach. Also note that the method only gives a lower-bound to the output's energy.
It *does not* give a guarantee that the output will have that energy, it can only estimate to which extent noise will deteriorate the quality of the outputs of the device.

The following flowchart explains the workings a bit better:



How does KnowYourLimits obtain the lower bound?
-----------------------------------------------

To fully understand the answer to that question, one should of course consult the original paper the method is based on. But roughly speaking, the method hinges
on computing certain entropic quantities related to the circui and the partition function of the Ising Hamiltonian. The mehod used for the first part (computing entropic quantities) is
computationally inexpensive given the compiled and routed circuit. On the other hand, computing the partition function is computationally more expensive. 

KnowYourLimits resorts to two algorithms to compute the partition functions: Markov Chain Monte Carlo and Tensor Networks. The user can then specify which method 
they prefer to compute the partition function. We shortly discuss their advantages and disadvantages below and refer to XXX for the different routines implemented.


Markov Chain Monte Carlo
*************************

Markov Chain Monte Carlo methods have the advantage that they are highly scalable. They can be used to obtain estimates for devices comprised of hundreds or thousands
of qubits within a resonable computational time. See Example XXX for a discussion of one case. However, they can only used efficiently and reliably to estimate partition functions
at relatively high temperatures. This means that the bounds obtained by KnowYourLlimits with this method are looser.

Tensor Networks
***************

Tensor network methods for computing partition functions have the advantage of reliably estimating the partition function for systems with up to a hundred qubits 
even at high temperatures. This means that they can be used to obtain tighter lower bounds on the energy of the circuit when compared to the Markov Chain methods.
However, it quickly becomes too expensive to use these methods for more than a hundred qubits or so.

as an input a noise model of an underlying device specified by a quantum channel. 
It is possible to specify noise models for both 1-qubit gates and 2-qubit gates, as many current implementations have significantly different error rates for these
two classes of gates. The quantum channels can be specified by their Kraus operators or a `Cirq noise model <https://quantumai.google/cirq/noise>`_, but we assume that
the noise is independent of the gate being implemeneted.


.. automodule:: main
    :members:


..  image:: ./images/know_limits_vs_google.png
    :align: center
    :width: 500
    :alt: Know your limits vs. Google