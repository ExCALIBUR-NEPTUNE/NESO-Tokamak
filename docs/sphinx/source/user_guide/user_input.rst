User Input
=========================

.. role:: xml(code)
   :language: xml

The details of each simulation are specified in two xml files: the mesh and the config.

Solver Specification
====================


Magnetic Field
---------------------

The magnetic field is specified as a function.

There are two options for the function name: "MagneticField" and "MagneticMeanField".

"MagneticField" indicates that the given components are evaluated as written.

"MagneticMeanField" indicates that the given components are evaluated as written on the x-y plane, and then rotated around the y-axis to generate an azimuthally symmetric 3D field.
For tokamaks, the poloidal plane is equated with the x-y plane.

.. code-block:: xml

    <FUNCTION NAME="MagneticField">
        <E VAR="Bx" VALUE="(y-yc)/sqrt((x-rc)*(x-rc)+(y-yc)*(y-yc))" />
        <E VAR="By" VALUE="-(x-rc)/sqrt((x-rc)*(x-rc)+(y-yc)*(y-yc))" />
        <E VAR="Bz" VALUE="0" />
    </FUNCTION>

Note that the :xml:`<E>` tag indicates that the quantity is evaluated from a symbolic expression.
This feature is also used within other tags.
Most simple mathematical operations and algebra can be used within such an expression, as well as previously specified parameters.
The arguments to the expression are the Cartesian coordinates x,y,z and simulation time t.



Fluid Specificiation
====================

Fluid species are listed under the :xml:`<SPECIES>` tag.
:xml:`<I>` and :xml:`<N>` are used for charged and neutral fluids respectively, and the :xml:`NAME` field must be filled.
For ions, mass and charge parameters must be specified, while for neutrals, only the mass must be specified.

The :xml:`FIELDS` field must also be filled with the fields that should be evolved for that species.
For instance: :xml:`<I NAME="D+" FIELDS="n,v,e">` creates a plasma species called "D+" with density, momentum and internal energy evolution enabled.

The specified fields must be initialised using the :xml:`<INITIAL>` tag:

.. code-block:: xml

    <I NAME="D+">
        <P> Mass = 2.0 </P>
        <P> Charge = 1.0 </P>
        <INITIAL>
            <E VAR="n" VALUE="sin(4*x)*sin(4*x)"/>
        </INITIAL>
    </KI>






Particle Specification
======================

PENKNIFE can use the NESO-Particles and NESO frameworks to allow users to use computational particles which interact with the fluid solver.

A Particle system is enabled using the :xml:`<PARTICLES>` field in the XML config file.

Global parameters for the particle system are specified in :xml:`<PARAMETERS>`

Particle information is specified on a per-species basis, within :xml:`<SPECIES>`

:xml:`<KN>` and :xml:`<KI>` are used for neutral particles and charged particles respectively, and the :xml:`NAME` field must be filled

Neutral species must have a mass parameter, and ions species must have both mass and charge parameters.
All species must also have an initial distribution specified, and may also optionally have sources and sinks.

.. code-block:: xml

    <KN NAME="D">
        <P> Mass = 2.0 </P>
        <INITIAL N="1000">
            <E VAR="n" VALUE="sin(4*x)*sin(4*x)"/>
        </INITIAL>
        <SOURCES>
            <P N="100">
                ...
            </P>
        </SOURCES>
        <SINK>
            ...
        </SINK>
        <BOUNDARYINTERACTION>
            ...
        </BOUNDARYINTERACTION>
    </KN>

Initial Conditions
------------------

The species initial conditions are specified in the :xml:`<INITIAL N="num">` field.

The density field can be specified by an analytic function for a probability distribution (everywhere between 0 and 1).
During initialisation, `num` particles are created using adaptive rejection sampling over this distribution.

.. code-block:: xml

    <INITIAL N="1000">
        <E VAR="W" VALUE="0.000005"/>
        <E VAR="n" VALUE="exp((-(x-rc-rs)*(x-rc-rs)/(sr*sr))+(-(y-yc-ys)*(y-yc-ys)/(sy*sy)))" />
        <E VAR="V" VALUE="0.001"/>
    </INITIAL>

Sources
-------

Sources are specified under a :xml:`<SOURCE N="num">` field.  Multiple sources can be created for each species.
Here, `num` indicates the number of particles added at each timestep.

The following xml indicates that there is a point source of particles at (2.75, 0), adding 100 particles of weight 0.05 at each time step, with thermal velocities at temperature 1eV

.. code-block:: xml

    <SOURCES>
        <P N="100">
            <E VAR="W" VALUE="0.05"/>
            <E VAR="X" VALUE="2.75"/>
            <E VAR="Y" VALUE="0"/>
            <E VAR="T" VALUE="1"/>
        </P>
    </SOURCES>

Sinks
-----

Sinks are specified under a :xml:`<SINK>` field.  Multiple sinks can be created for each species.
The specified function indicates the probability that a particle is removed at each timestep.
The probability function is evaluated at each particle's position, which is then marked for removal with that probability.

The following xml indicates a sink with probability with a moving Gaussian.

.. code-block:: xml

    <SINK>
        <E VAR="n" VALUE="exp((-(x-2-0.75*sin(150*t))*(x-2-0.75*sin(150*t))/(0.01))+(-(y+0.75*cos(150*t))*(y+0.75*cos(150*t))/(0.01)))"/>
    </SINK>


Boundary Conditions
-------------------

The particles may interact with the boundary regions by either reflecting, being absorbed, or taking part in a surface reaction.

Optionally, the user can use the Reactions/VANTAGE framework by populating the :xml:`<VANTAGE>` XML field.

The three currently supported reactions are ionisation, recombination and charge exchange.


Reactions coupling
==================

Both bulk and surface reactions can be added under the :xml:`<VANTAGE>` XML field:
:xml:`<R>` indicates that the reaction applied to all particles.
:xml:`<S>` indicates that the reaction is only applied to particles that hit a boundary.

Bulk Reactions
--------------

Ionisation
~~~~~~~~~~

.. code-block:: xml

    <R TYPE="Ionisation" SPECIES="D">
        <RATE TYPE="AMJUEL" VALUE="1"/>
        <CROSSSECTION TYPE="Constant" VALUE="1.0"/>
    </R>

Recombination
~~~~~~~~~~~~~

.. code-block:: xml

    <R TYPE="Recombination" SPECIES="D">
        <RATE TYPE="AMJUEL" VALUE="1"/>
        <CROSSSECTION TYPE="Constant" VALUE="1.0"/>
    </R>

Charge Exchange
~~~~~~~~~~~~~~~

.. code-block:: xml

    <R TYPE="ChargeExchange" SPECIES="D">
        <RATE TYPE="AMJUEL" VALUE="1"/>
        <CROSSSECTION TYPE="Constant" VALUE="1.0"/>
    </R>

Surface Reactions
-----------------

For each surface reaction the species it involves are listed after the type.
It is also possible to have separate xml fields for each species and type of reaction if different parameters are desired.


The :xml:`<REGION>` field has references to the IDs of the boundary regions at which the boundary reaction should apply.
The boundary regions are the same ones identified in the :xml:`<BOUNDARYREGIONS>` field inside :xml:`<CONDITIONS>`.
Each boundary region has a surcface boundary controller which handles all the reactions that occur on that region.

Specular Reflection
~~~~~~~~~~~~~~~~~~~

The :xml:`VALUE` field indicates the relative proportion of weight that is carried by the reflected particle.


.. code-block:: xml
    
    <S TYPE="Specular" SPECIES="D,T">
        <REGION REF="1,2"/>
        <RATE TYPE="Constant" VALUE="0.9"/>
    </S>


Absorption
~~~~~~~~~~

The :xml:`VALUE` field indicates the relative proportion of weight that is absorbed.

.. code-block:: xml
    
    <S TYPE="Absorption" SPECIES="D,T">
        <REGION REF="1,2"/>
        <RATE TYPE="Constant" VALUE="0.1"/>
    </S>
    

Thermal Reflection
~~~~~~~~~~~~~~~~~~

The :xml:`VALUE` field indicates the relative proportion of weight that is carried by the reflected particle.


.. code-block:: xml

    <S TYPE="Thermal" SPECIES="D,T">
        <REGION REF="1,2"/>
        <RATE TYPE="Constant" VALUE="0.5"/>
    </S>