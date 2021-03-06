#############################
Partition function estimators
#############################
The main methods to compute partition functions. To see the details of how the partition function is computed with tensor network methods, we refer to the `quimb documentation <https://quimb.readthedocs.io/en/latest/_autosummary/quimb.tensor.tensor_gen.html#quimb.tensor.tensor_gen.HTN_classical_partition_function_from_edges>`_.
But the important message is that they can compute the partition function relatively efficiently for system sizes up to 70 qubits as long as the underlying graph is sufficiently sparse.


For the Markov chain method, a more detailed explanation may be found in the corresponding page in the auxilliary modules, see :ref:`Markov chain methods`. Their main advantage is that they can compute partition functions for significantly larger system sizes, up to thousands of qubits. 
However, they generally obtain loser estimates and for smaller sizes the tensor networks are faster.

.. automodule:: partfunc_estimators
    :members: