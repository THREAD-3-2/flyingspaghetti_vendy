.. _formulation:


===============================
Velocity based beam formulation
===============================

The source code `flyingspaghetti_vendy <https://github.com/THREAD-3-2/flyingspaghetti_vendy>`_  is based on novel energy preserving velocity based formulation for the geometrically exact beams  proposed by `Zupan and Zupan, (2019) <https://doi.org/10.1007/s11071-018-4634-y>`_.
For the sake of brevity, the final discretized governing equations are presented here. The time discretization employed here is based on
midpoint rule while the spatial discretization is based on Galerkin finite element method. The final discretized momentum balance equations are:

.. math::
       :name: eq:1

       \begin{aligned}
                        \int_{0}^{L} \left[\frac{\rho A}{h}\left(\boldsymbol{v}^{\left[n+1\right]} - \boldsymbol{v}^{\left[n\right]}\right) P_{i} + \boldsymbol{n}^{\left[n+1/2\right]}P_{i}^{\prime} - \boldsymbol{\tilde{n}}^{\left[n+1/2\right]}P_{i}\right] \,\mathrm{dx} - \delta_{p}\boldsymbol{f}_{e}^{\left[n+1/2\right]} = \boldsymbol{0}
       \end{aligned}

.. math::
       :name: eq:2

       \begin{aligned}
                    &\int_{0}^{L} \bigg[\frac{\boldsymbol{J}_{\rho} }{h}\left(\boldsymbol{\mathit{\Omega}}^{\left[n+1\right]} - \boldsymbol{\mathit{\Omega}}^{\left[    n\right]}\right) P_{i} + \boldsymbol{\overline{\mathit{\Omega}}}\times \boldsymbol{J}_{\rho}\boldsymbol{\overline{\mathit{\Omega}}}P_{i} - \boldsymbol{K}^{\left[n+1/2\right]} \times\boldsymbol{M}^{\left[n+1/2\right]}P_{i} +\boldsymbol{M}^{\left[n+1/2\right]}P_{i}^{\prime}\\& -  \left(\boldsymbol{\mathit{\Gamma}}^{\left[n+1/2\right]} - \boldsymbol{\mathit{\Gamma}}_{\mathit{0}}\right) \times \boldsymbol{N}^{\left[n+1/2\right]} P_{i}  -\left(\boldsymbol{\hat{q}}^{\ast\left[n+1/2\right]} \circ \boldsymbol{\tilde{m}}^{\left[n+1/2\right]} \circ \boldsymbol{\hat{q}}^{\left[n+1/2\right]}\right)P_{i}\bigg] \,\mathrm{dx} -\delta_{p}\boldsymbol{M}_{e}^{\left[n+1/2\right]} = \boldsymbol{0}
       \end{aligned}

where :math:`\boldsymbol{v},\boldsymbol{\Omega}` are the velocities and angular velocities in fixed and local basis respectively, :math:`\rho` and :math:`\boldsymbol{J}_{\rho}` are material density and mass moment of inertia, :math:`\boldsymbol{N}` and :math:`\boldsymbol{M}` are respectively the stress resultant forces and moments, :math:`\boldsymbol{\Gamma}` and :math:`\boldsymbol{K}` are the vector of strain and curvature and :math:`\boldsymbol{\tilde{n}}` and :math:`\boldsymbol{\tilde{m}}` are the vector of distributed force and distributed moments. In the above equations, the quantities in the fixed basis are denoted in lower case while the ones in local basis are denoted with upper case notations.

In the present formulation, rotational quaternions are employed for the description of three-dimensional rotations. However, in contrast to many existing approaches, the interpolated variables here are the velocities and angular velocities at midtime, :math:`(\boldsymbol{\overline{v}},\boldsymbol{\overline{\Omega}})`. The interpolation points :math:`x_i (i= 1,2,3,\dots,p)` are chosen equidistantly from the interval :math:`[0,L]` with :math:`x_0=0` and :math:`x_p = L`. The crucial idea here is that the angular velocities are locally additive quantities and hence the angular momentum balance equation is expressed here in local basis. Thus standard additive type interpolation functions can be used while avoiding the multiplicative update or special transformations for the rotational increments.

.. math::
       :name: eq:3

       \begin{aligned}
            \boldsymbol{\overline{v}}(x,t) &= \sum_{i=1}^{p} P_{i}(x)\boldsymbol{\overline{v}}^{i}(t)
       \end{aligned}

.. math::
       :name: eq:4

       \begin{aligned}
            \boldsymbol{\overline{\mathit{\Omega}}}(x,t) &= \sum_{i=1}^{p} P_{i}(x)\boldsymbol{\overline{\mathit{\Omega}}}^{i}(t)
       \end{aligned}

The computational benefits of the approach here are: (i) the additivity of the primary unknowns enables the use of standard interpolation techniques in space and simplifies the iteration procedure; (ii) the energy of the conservative systems is preserved without the need to introduce any special approxmiation for the roation or strain field and (iii) the kinematic compatibility equations are satisfied with the same level of accuracy as the governing equations. For more detailed understanding, the user is suggested to go through the paper.
