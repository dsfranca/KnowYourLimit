###############################
Simple estimate generator
###############################

In this module we provide some functions that make it easy to perform some simple estimates "out of the box". So far it mostly focus on QAOA circuits and quantum annealing algorithms.

We provide functions that, given an architecture, noise model and Hamiltonian we wish to minimize, allow for lower-bounding the energy of the output of circuits or annealers. 

We provice functions to lower-bound the energy on SK models and random d-regular graphs for square architectures under depolarizing noise. These allow for an easy testing of how limited connectivity of quantum devices affects the performance of the QAOA under noise.



.. automodule:: simplified_estimator
    :members: