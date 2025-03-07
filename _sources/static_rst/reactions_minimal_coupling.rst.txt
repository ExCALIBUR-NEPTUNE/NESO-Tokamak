=================================================
NESO-Tokamak-Reactions coupling minimal equations
=================================================

Single field anisotropic diffusion with sources and outflow
-----------------------------------------------------------

This can be in 2D or 3D.

.. math:: \frac{\partial n}{\partial t} + \nabla \cdot \vec{\Gamma} = S_n

with the diffusive flux given by :math:`\vec{\Gamma} = D \cdot \nabla n`
with D being the anisotropic diffusion tensor
:math:`D = \vec{b}\vec{b} k_\parallel + (I - \vec{b}\vec{b}) k_\perp`
where the :math:`k_\parallel` and :math:`k_\perp` should have the option
of being functions of fields, and :math:`\vec{b}` is the local unit
vector tangential to the magnetic field.

The source :math:`S_n` should be given by Reactions, with the idea of
starting off with ionisation only, so it would always be positive in
this example. For the purpose of Reactions and the Bohm BC, we will
assume an isothermal plasma here.

The boundary condition should be set to Neumann in the directions
perpendicular to the magnetic field, while the parallel component of the
flux at the boundary

.. math:: \vec{\Gamma}\cdot \vec{b} = n c_s

where :math:`c_s` in this case is a constant (isothermal Bohm speed
:math:`c_s = \sqrt{kT/m_i}`).

Once basic coupling has been done with this, we can add another field.

Two diffusive fields with particle and energy coupling with Reactions
---------------------------------------------------------------------

Same as above for :math:`n`-field, with the addition of non-constant
:math:`T`.

We will ignore all advective and field terms for the first pass, and
write a pressure equation as

.. math:: \frac{3}{2} \frac{\partial p}{\partial t} + \nabla \cdot \vec{q} = S_E + Q,

where :math:`p = nkT` and

.. math:: \vec{q} = \kappa \cdot \nabla T

with the conductivity tensor given as

.. math:: \kappa =  \kappa_\parallel(T) \vec{b}\vec{b} + \kappa_\perp(T) (I-\vec{b}\vec{b})

with parallel and perpendicular conductivities being functions of
temperature to start off with.

Energy sinks :math:`S_E` should again be supplied by Reactions, with
some heating :math:`Q` in the core being an input parameter.

At the parallel boundaries, we set the energy outflow to

.. math:: \vec{q}\cdot \vec{b} = \gamma nkTc_s,

where now :math:`c_s=\sqrt{kT/m_i}` for both this and the particle BC,
and :math:`\gamma \approx 7`. No outflow/Neumann at perpendicular
boundaries.

Only after both of these models have been tested with Reactions coupling
does it make sense to move towards a model with an explicit momentum
equation, such as a reduced version of Hermes-3 mean-field equations.
