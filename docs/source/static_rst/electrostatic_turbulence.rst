Multi-species electrostatic turbulence system with neutrals
-----------------------------------------------------------
Solves a system of :math:`4+3S` fields representing :math:`S` ion species.
Electrons are always present, and at least one ion species must be present.
Explicitly the fields are:

.. math:: \left[\omega, n_e,  m_e n_e v_{\parallel e}, \frac{3}{2}p_e, n_i m_i n_i v_{\par i}, \frac{3}{2}p_i, ...]

The continuity equations are:

.. math:: \frac{\partial n}{\partial t} = \del \cdot \left[n\left(\mathbf{v}_{\mathbf{E}\times\mathbf{B}} + \mathbf{b}v_{\parallel e} + \mathbf{v}_de\right)\right] + S_n
.. math:: \frac{\partial \omega}{\partial t}

.. math:: \frac{\partial}{\partial t} \left(m_e n_e v_{\parallel e}\right) = \del \cdot \left[m_e n_e v_{\parallel e}\left(\mathbf{v}_{\mathbf{E}\times\mathbf{B}} + \mathbf{b}v_{\parallel e} + \mathbf{v}_de\right)\right] - \mathbf{b}\cdot\del p_e - enE_{\parallel} - F_{ei}
.. math:: \frac{\partial}{\partial t} \left(\frac{3}{2}p_e\right) = \del \cdot \left[\frac{3}{2}p_e\left(\mathbf{v}_{\mathbf{E}\times\mathbf{B}} + \mathbf{b}v_{\parallel e} + \frac{5}{2}\mathbf{v}_de\right)\right] - p_e\del\cdot\left(\mathbf{v}_{\mathbf{E}\times\mathbf{B}} + \mathbf{b}v_{\parallel e} right) + \del\cdot\left(\kappa_{\parallel e}\mathbf{b}\mathbf{b}\cdot\del T_e\right) + S_{Ee} + W_{ei}
  
.. math:: \frac{\partial}{\partial t} \left(m_i n_i v_{\parallel i}\right) = \del \cdot \left[m_i n_i v_{\parallel i}\left(\mathbf{v}_{\mathbf{E}\times\mathbf{B}} + \mathbf{b}v_{\parallel i} + \mathbf{v}_di\right)\right] - \mathbf{b}\cdot\del p_i + ZenE_{\parallel} + F_{ei}
.. math:: \frac{\partial}{\partial t} \left(\frac{3}{2}p_i\right) = \del \cdot \left[\frac{3}{2}p_i\left(\mathbf{v}_{\mathbf{E}\times\mathbf{B}} + \mathbf{b}v_{\parallel i} + \frac{5}{2}\mathbf{v}_di\right)\right] - p_i\del\cdot\left(\mathbf{v}_{\mathbf{E}\times\mathbf{B}} + \mathbf{b}v_{\parallel i} right) + \del\cdot\left(\kappa_{\parallel i}\mathbf{b}\mathbf{b}\cdot\del T_i\right) + S_{Ei} + S_n\frac{1}{2}m_i n_i v_{\parallel i}^2 - W_{ei} + \frac{p_i}{en_0}\del\cdot\left(\mathbf{J}_{\parallel}+\mathbf{J}_d\right)

The last two equations are duplicated for each ion species.


Finally the vorticity equation is:
.. math:: \frac{\partial\omega}{\partial t} = \del \cdot \left[ \right]

Boundary Conditions
-------------------

The boundary conditions for this system of equations are as follows.

For the electron density

.. math:: \frac{\partial n}{\partial \mathbf{n}} = \mp \frac{n}{c_s} \frac{\partial v_{\parallel i}{\partial \mathbf{n}}

For parallel electron velocity

.. math:: v_{\parallel e} = \pm c_s \exp\left(\Lambda - \frac{\phi}{\Te}\right)

For the electron energy

.. math:: \frac{\partial \frac{3}{2}p_e}{\partial \mathbf{n}} = \gamma_e n_e T_e c_s

For parallel ion velocity

.. math:: v_{\parallel i} = \pm c_s

For the ion energy

.. math:: \frac{\partial \frac{3}{2}p_i}{\partial \mathbf{n}} = \gamma_i n_i T_i c_s

For vorticity

.. math:: \omega = -\cos^2\alpha\left(\frac{\partial v_{\parallel i}{\partial \mathbf{n}}\right)^2 \pm c_s \frac{\partial^2 v_{\parallel i}{\partial \mathbf{n}^2}