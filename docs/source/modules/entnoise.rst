##############################
Noise and entropic convergence
##############################

.. _my-reference-label:

Quantum channels
------------------------
One of the basic objects required to apply our technique is a quantum channel modelling the noise affecting the quantum device. We will depart from the assumption that the noise affecting the device is identical and independently distributed on each qubit and that we have one qubit noise.

We will later expand our package to also allow for each qubit having a different quantum channel affecting it and also to allow for different noise rates for 1-qubit and 2-qubit gates.

We will use `Cirq quantum channels <https://quantumai.google/cirq/noise>`_ to specify the noise models. 

We automatically compute the fixed point of the channel and the contraction coefficient of the relative entropy with respect to the fixed point. This is the main information on the channel we require to compute our lower bound.

Recall that the technique only works if the quantum channel has a unique fixed point and it is not pure. You will obtain warnings if the quantum channel specified does not satisfy these conditions.

Lindbladians
~~~~~~~~~~~~~
Quantum annealers being computational devices with continuous time, it is also more natural to model the noise affecting them in continuous time as well.
That is, we will model the noise with a Lindbladian.

Thus, KnowYourLimits also has functions to specify a noise model in continuous time. 
In the current implementation, it is assumed that the quantum annealer is implementing a time-dependent evolutions, but that the noise model itself is time-independent.
Furthermore, we will assume that the Lindbladian affects each qubit independently and that the noise affecting all qubits is the same.
For a more detailed discussion of this noise model, we refer to Sec. ID of the supplementary material of `"Limitations of optimization algorithms on noisy quantum devices" <https://doi.org/10.1038/s41567-021-01356-3>`_.

Unlike it was the case for discrete time quantum channels, Cirq does not readily allow the specification of a Lindbladian :math:`\mathcal{L}`.
However, many Lindlbadians of physical interest are of the form :math:`\mathcal{L}=r(T-\textrm{id})`, for some quantum channel :math:`T` and a noise rate :math:`r>0`. 
Thus, the Lindbladian class will accept specifying both the Lindbladian in terms of a quantum channel and a rate or in terms of the process matrix of the Lindbladian.

.. automodule:: noise_model
    :members:


Entropic convergence for circuits
-----------------------------------
This module contains the function used to estimate the relative entropy between the output of a circuit and the fixed point of the noise affecting the circuit.

This estimate is performed using the entropic convergence result stated in Lemma 3 of `"Limitations of optimization algorithms on noisy quantum devices" <https://doi.org/10.1038/s41567-021-01356-3>`_. As discussed in more detailed in the aforementioned paper, this estimate is particularly efrfective for QAOA-like circuits or whenever the fixed point of the noise is close to maximally mixed.

.. automodule:: circuit_analyser
    :members:



Entropic convergence for annealers
-----------------------------------
This module contains the functions used to estimate the relative entropy between the output of a noisy quantum annealer and the fixed point of the noise.

This estimate is performed using the entropic convergence result stated in Theorem 1 of `"Limitations of optimization algorithms on noisy quantum devices" <https://doi.org/10.1038/s41567-021-01356-3>`_.

Note that we assume that the underlying annealing path is of the form :math:`f(t/T)H_0+g(t/T)H_I`, where :math:`H_0=-\sum_i X_i`. This reflects the annealing paths that can be implmeented in contemporary platforms.

.. automodule:: annealer_analyser
    :members:
