======================
Particle Specification
======================

NESO-Tokamak can use the NESO-Particles and NESO frameworks to allow users to use computational particles which interact with the fluid solver.

A Particle system is enabled using the `<PARTICLES>`field in the XML config file.

Global parameters for the particle system are specified in `<PARAMETERS>`

Particle information is specified on a per-species basis, within `<SPECIES>`

Each species must have mass and charge fields, as well as initial conditions, and may optionally have sources and sinks.

Initial Conditions
------------------

The species initial conditions are specified in the `<INITIAL N="num">`field.

The density field can be specified by an analytic function for a probability distribution.
During initialisation, `num` particles are created using adaptive rejection sampling over this distribution.

Sources
-------

Sources are specified under a `<SOURCE N="num">`field.  Multiple sources can be created for each species.
Here, `num` indicates the number of particles added at each timestep.

Sinks
-----

Sinks are specified under a `SINK` field.  Multiple sinks can be created for each species.
The specified function indicates the probability that a particle is removed at each timestep.
The probability function is evaluated at each particle's position, which is then marked for removal with that probability.



