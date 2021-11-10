#############################
Markov chain methods
#############################

This module implements functions required to estimate the partition function of an Ising model with Markov Chain Monte Carlo methods. . 

To compute the partition function, we resort to the telescopic product:

.. math::

   Z_{\beta_{k}}=Z_{\beta_{0}}\frac{Z_{\beta_{1}}}{Z_{\beta_{0}}}\frac{Z_{\beta_{2}}}{Z_{\beta_{1}}}\cdots\frac{Z_{\beta_{k}}}{Z_{\beta_{k-1}}}

Where :math:`\beta_0=0` and the choice of the :math:`\beta_i` is the annealing schedule. To estimate each ratio we can resort to the fact that the expectation value of the random variable 

.. math::

   \textrm{exp}(-(\beta_{i+1}-\beta_{i})H(x))

with respect to the Gibbs state :math:`e^{-\beta_{i}H}/Z_{\beta_{i}}` is given by the ratio :math:`Z_{\beta_{i+1}}/Z_{\beta_{i}}`. By picking the :math:`\beta_{i}` close enough to each other, we can reliably estimate the ration by sampling from that Gibbs state and computing the empirical mean.

A huge amount of literature is dedicated to how to sample from the underlying Gibbs states as efficiently as possible. The same can be said about how to pick the annealing strategy so that the number of different :math:`\beta_{i}` is as small as possible.

The first version of know your limits only implements the most basic version of these two subroutines (sampling from the Gibbs state and the annealing schedule).

For the sampling we use the Metropolis/Glauber dynamics algorithm. And for the annealing we a uniform schedule. That is, the :math:`\beta_{i}` are always incremented by the same amount. However, the user can also specify the annealing schedule, so the code can easily handle more sophisticated schedules the user wants to implement. In future versions we plan to implement more advanced sampling routines and annealing schedules.

Also note that the Markov chain method is only (provably) efficient for :math:`\beta_{k}=\mathcal{O}(\|A\|^{-1})`. For that range, a number of Metropolis steps of order :math:`\mathcal{O}(n\log(n))` should suffice to obtain a sample that is close to the Gibbs state. See e.g. the book `Markov Chains and Mixing Times <https://pages.uoregon.edu/dlevin/MARKOV/markovmixing.pdf>`_ for a thorough discussion of all the topics discussed above.


.. automodule:: mcmc
    :members: