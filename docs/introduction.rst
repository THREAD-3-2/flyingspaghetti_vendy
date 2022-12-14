.. _introduction:

=========
Example 
=========

Problem description 
==========================

This is a classical example introduced by `Simo and Vu-Quoc, (1988) <https://doi.org/10.1016/0045-7825(88)90073-4>`_ and it involves extremely large displacements and rotations and served later also for the demonstration of numerical stability of energy preserving algorithms. The beam lies initially inclined in :math:`XZ`-plane and is subjected at the lower end to point force :math:`f_{x}` and point moments :math:`h_{y}` and :math:`h_{z}` with constant direction and piecewise linear magnitude in time as shown in the figure.

.. figure:: /figures/flying_spaghetti.png

   Figure 1. Free flight of a flexible beam

At :math:`t = 5`, the point loads vanish leaving the beam to fly in
free motion until the chosen final time. The remaining parameters of the beam
are: :math:`EA = GA_{2} = GA_{3} = 10000`,
:math:`GI_{1} = EI_{2} = EI_{3} = 500`,
:math:`\rho A = 1, \boldsymbol{J}_{\rho} = \text{diag}[10, 10, 10]`.
