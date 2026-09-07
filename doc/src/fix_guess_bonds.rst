.. index:: fix guess_bonds

fix guess_bonds command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID guess_bonds Nevery Nrepeat Nfreq bond_order_cutoff mode values

* ID, group-ID are documented in :doc:`fix <fix>` command
* guess_bonds = style name of this fix command
* Nevery = accumulate atom attributes once every this many steps
* Nrepeat = # of times to accumulate atom attributes
* Nfreq = make stored atom attributes (history) available every this many steps
* fraction = add bond if a bond was present for at least this fraction of sampled timesteps (value between 0.0 and 1.0)
* mode can be either *distance* or *pauling*

  .. parsed-literal::

       *distance* values = prefactor args
         prefactor = tolerance factor when guessing bonds
         args = list of atom type/radius pairs
       *pauling* values = BOcut args
         BOcut = bond order cutoff to decide if bond exists
         args = list of each atom type pair followed by their equilibrium bond distance and bond softness parameter

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all guess_bonds 10 10 100 0.5 distance 1.15 c1 0.70 hc 0.37 oc 0.6 n 0.65
   fix 1 all guess_bonds 10 10 100 0.5 &
     pauling 0.3 C C 1.54 0.353 C H 1.09 0.399 C O 1.43 0.279 C N 1.47 0.298 O N 1.44 0.26

Description
"""""""""""

Guess if a bond exists between each pair of atoms in the simulation. If
a bond is determined to exist, it is added to the internal bond list.
This fix is meant to be used with atom styles that do not use bond lists
when calculating the physics of the simulation, such as reactive force
fields or machine-learned interatomic potentials. However, a bond style
must be defined to use this fix. The :doc:`bond_style zero <bond_zero>`
can be used to prevent the added bonds from affecting the dynamics.
Other fixes or computes that used the bond list or act on molecules,
such as :doc:`compute fragment/atom <compute_cluster_atom>`,
:doc:`compute chunk/atom molecule <compute_chunk_atom>`, or :doc:`fix
evaporate <fix_evaporate>` can then be applied for analysis or
algorithmic purposes.

The *Nevery*, *Nrepeat*, and *Nfreq* values are used for an internally-
created :doc:`fix store/state <fix_store_state>` command that saves bond
history. If the fraction of times that is a bond is identified, out of
*Nrepeat* sampled timesteps, is greater than *fraction*, then this fix
assigns a bond. For example if *fraction* = 0.5, and a bond was
identified between atoms *i* and *j* for 60 of the 100 previously
sampled timesteps, then a bond is added between atom *i* and *j* on the
current timestep. The bonds identified by this fix are added as explicit
bonds so that the geometry of bond pairs is accessible to other
commands. A bond type of 1 is assigned to all bonds added by this fix.

.. note::

   This fix should generally not be used with models that require
   explicit bonds.  This fix adds and *deletes* bonds but does not
   update angles, dihedrals, or impropers.  This fix is meant to assist
   with on-the-fly analysis for atom styles that do not define explicit
   bonds, such as machine learning interatomic potentials or bond-order
   force fields. :doc:`Bond_style zero <bond_zero>` can be used to
   prevent these guessed bonds from affecting the dynamics of the
   simulation.

Bond order is calculated differently depending on whether the *distance*
or *pauling* mode is specified.

In *distance* mode, if a bond satisfies the distance cutoff criterion on
a given timestep, it is assigned a bond order of 1.0, otherwise it is
assigned a bond order of 0.0. In *distance* mode, the distance cutoff is
defined by the following equation:

.. math::

   cutoff = prefactor*(r_1 + r_2)

where :math:`r_1` is the atomic radius of the first atom and :math:`r_2`
is the atomic radius of the second atom, and a bond order of 1.0 is
assigned between these two atoms if their interatomic distance is less
than the cutoff. All existing bonds that do not fulfill this criterion
are removed. The prefactor must be provided directly after the
*distance* keyword. Atomic radii must be provided for every atom type
following the prefactor, where the numeric atom type or type label is
followed by the atomic radius (see example above).

In *pauling* mode, the bond order is calculated with the following equation:

.. math::

   BO = \exp\left(\frac{r_0 - r}{b}\right)

where BO is bond order, :math:`r_0` is the equilibrium bond length of a
single bond (as opposed to a double bond or triple bond) between the two
atoms, :math:`r` is the instantaneous distance between two atoms, and
:math:`b` is a bond softness parameter. Values for :math:`r_0` and
:math:`b` must be provided for every pair of atom types, where the
numeric atom types or type labels are followed by :math:`r_0` and
:math:`b` (see example above). An explicit bond is added between these
two atoms if their calculated bond order is less than
*bond_order_cutoff*. All existing bonds that do not fulfill this
criterion are removed.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Bonds are always written to restart files.  No additional information
about this fix is written to :doc:`binary restart files <restart>`.
None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
none

Related commands
""""""""""""""""

:doc:`bond_style zero <bond_zero>`, :doc:`compute mliap <compute_mliap>`

Default
"""""""

none
