###############################
Simple estimate generator
###############################

In this module we provide some functions that make it easy to perform some simple estimates "out of the box". So far it mostly focus on QAOA circuits.

In particular, we provice functions to lower-bound the energy on SK models and random d-regular graphs for square architectures under depolarizing noise. These allow for an easy testing of how limited connectivity of quantum devices affects the performance of the QAOA under noise.

.. automodule:: simplified_estimator
    :members: