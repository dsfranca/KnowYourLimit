###############################
Compilation and routing
###############################

Whenever implementing a quantum circuit on an given quantum device, there are two limitations of any real device that have to be taken into account:
 
 #. The device cannot implement arbitrary 1-qubit and 2-qubit gates. Thus, the gates in the circuit have to be compiled into a circuit that only containts native gates.
 #. The device has a limited connectivity. That is, we cannot implement gates between arbitrary 2-qubits, but only those specified by some connectivity graph.
 
Both of these points increase the physical depth of a circuit when implementing it on near-term devices. Although both of these points only incur in a polynomial increase of the depth of the circuit and, thus, is not a major issue for fault-tolerant quantum computation, they can make a significant difference in the absence of error correction.

The functions in this module allow for both compling a given circuit into a native gate set and to perform the routing of the qubits in efficient ways. To perform the compilation we will resort to `Cirq <https://quantumai.google/cirq/circuits>`_ and to perform the routing we will resort to `TKET <https://cambridgequantum.com/tket/>`_ (see also `this tutorial <https://quantumai.google/cirq/experiments/qaoa/routing_with_tket>`_ for a discussion on how to use tket with cirq).

This first version will mostly focus on QAOA circuits and we will use some of the functions developed in the `Recirq package.<https://github.com/quantumlib/ReCirq>`_

.. automodule:: circuit_compiler
    :members: