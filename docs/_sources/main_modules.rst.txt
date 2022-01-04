Main Modules
=============

The package contains three main technical modules, outlined below.

Partition function estimators
-----------------------------
Functions required to estimate partition functions of Ising models with tensor networks or Markov chain Monte Carlo.


..  toctree::

    modules/partfunc


Compile and route
-----------------
Routines to compile and route quantum circuits to various architectures.

..  toctree::

    modules/comproute

Noise and entropic convergence
-------------------------------
Routines to estimate the entropic convergence of various noise models.

..  toctree::

    modules/entnoise



We also have one module that does simple estimates without all the bells and whistles of KnowYourLimits:

Simple estimates
----------------
Generates some simple estimates without having to estimate all the relevant parameters.

..  toctree::

    modules/simple_est