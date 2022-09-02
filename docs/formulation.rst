.. _formulation:


===============================
Velocity based beam formulation
===============================

The source code flyingspaghetti_vendy,
https://github.com/THREAD-3-2/flyingspaghetti_vendy.git\  is based on novel energy preserving velocity based formulation for the geometrically exact beams (Zupan,
Zupan, 2018) https://doi.org/10.1007/s11071-018-4634-y\. For the
sake of brevity, the final discretized governing equations are presented
here while for more detailed understanding, the user is suggested to go
through the paper. The time discretization employed here is based on
midpoint rule while the spatial descritization is based on Galerkin
finite element method.

.. math::
       \begin{aligned}
                        \int_{0}^{L} \left[\frac{\rho A}{h}\left(\boldsymbol{v}^{\left[n+1\right]} - \boldsymbol{v}^{\left[n\right]}\right) P_{i} + \boldsymbol{n}^{\left[n+1/2\right]}P_{i}^{\prime} - \boldsymbol{\tilde{n}}^{\left[n+1/2\right]}P_{i}\right] \,\mathrm{dx} - \delta_{p}\boldsymbol{f}_{e}^{\left[n+1/2\right]} = \boldsymbol{0}
       \end{aligned}

.. math::
       \begin{aligned}
                    &\int_{0}^{L} \bigg[\frac{\boldsymbol{J}_{\rho} }{h}\left(\boldsymbol{\mathit{\Omega}}^{\left[n+1\right]} - \boldsymbol{\mathit{\Omega}}^{\left[    n\right]}\right) P_{i} + \boldsymbol{\overline{\mathit{\Omega}}}\times \boldsymbol{J}_{\rho}\boldsymbol{\overline{\mathit{\Omega}}}P_{i} - \boldsymbol{K}^{\left[n+1/2\right]} \times\boldsymbol{M}^{\left[n+1/2\right]}P_{i} +\boldsymbol{M}^{\left[n+1/2\right]}P_{i}^{\prime}\\& -  \left(\boldsymbol{\mathit{\Gamma}}^{\left[n+1/2\right]} - \boldsymbol{\mathit{\Gamma}}_{\mathit{0}}\right) \times \boldsymbol{N}^{\left[n+1/2\right]} P_{i}  -\left(\boldsymbol{\hat{q}}^{\ast\left[n+1/2\right]} \circ \boldsymbol{\tilde{m}}^{\left[n+1/2\right]} \circ \boldsymbol{\hat{q}}^{\left[n+1/2\right]}\right)P_{i}\bigg] \,\mathrm{dx} -\delta_{p}\boldsymbol{M}_{e}^{\left[n+1/2\right]} = \boldsymbol{0}
       \end{aligned}

The interpolated variables are velocities and angular velocities at
midtime
:math:`(\boldsymbol{\overline{v}},\boldsymbol{\overline{\Omega}})`.
The interpolation points :math:`x_i; i= 1,2,3,\dots,p` are chosen
equidistantly from the interval :math:`[0,L]` with :math:`x_0=0` and
:math:`x_p = L`. The crucial idea here is that the angular velocities
are locally additive quantities and hence the angular momentum balance
equation is expressed here in local basis and thus standard additive
type interpolation functions can be used and avoids the multiplicative
update or special transformations for rotational increments.

.. math::
       \begin{aligned}
            \boldsymbol{\overline{v}}(x,t) &= \sum_{i=1}^{p} P_{i}(x)\boldsymbol{\overline{v}}^{i}(t)\\
            \boldsymbol{\overline{\mathit{\Omega}}}(x,t) &= \sum_{i=1}^{p} P_{i}(x)\boldsymbol{\overline{\mathit{\Omega}}}^{i}(t)
       \end{aligned}