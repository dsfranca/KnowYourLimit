##############################
Noise and entropic convergence
##############################


Quantum channels
------------------------
One of the basic objects required to apply our technique is a quantum channel modelling the noise affecting the quantum device. We will depart from the assumption that the noise affecting the device is identical and independently distributed on each qubit and that we have one qubit noise.

We will later expand our package to also allow for each qubit having a different quantum channel affecting it and also to allow for different noise rates for 1-qubit and 2-qubit gates.

We will use `Cirq quantum channels <https://quantumai.google/cirq/noise>`_ to specify the noise models. 

We automatically compute the fixed point of the channel and the contraction coefficient of the relative entropy with respect to the fixed point. This is the main information on the channel we require to compute our lower bound.

Recall that the technique only works if the quantum channel has a unique fixed point and it is not pure. You will obtain warnings if the quantum channel specified does not satisfy these conditions.

.. automodule:: noise_model
    :members:


Entropic convergence
------------------------
This module contains the function used to estimate the relative entropy between the output of a circuit and the fixed point of the noise affecting the circuit.

.. automodule:: circuit_analyser
    :members: