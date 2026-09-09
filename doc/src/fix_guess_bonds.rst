.. index:: fix guess_bonds

fix assign/bonds command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID assign/bonds Nevery Nrepeat Nfreq fraction c_ID

* ID, group-ID are documented in :doc:`fix <fix>` command
* assign/bonds = style name of this fix command
* Nevery = accumulate atom attributes once every this many steps
* Nrepeat = # of times to accumulate atom attributes
* Nfreq = make stored atom attributes (history) available every this many steps
* fraction = add explicit bond if bond was guessed for at least this fraction of sampled timesteps (value between 0.0 and 1.0)
* c_ID = 'c_' followed by the ID of a previously defined :doc:`compute guess/bonds <compute_guess_bonds>`

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all assign/bonds 10 10 100 0.5 c_guess_bonds

Description
"""""""""""

Add explicit bonds to the simulation if they were identified by
:doc:`compute guess/bonds <compute_guess_bonds>` and satisfy the
time-averaging criteria. If a bond is determined to exist, it is added to
the internal bond list. This fix is meant to be used with atom styles
that do not use bond lists when calculating the physics of the
simulation, such as reactive force fields or machine-learned interatomic
potentials. However, a bond style must be defined to use this fix. The
:doc:`bond_style zero <bond_zero>` can be used to prevent the added bonds
from affecting the dynamics. Other fixes or computes that used the bond
list or act on molecules, such as :doc:`compute fragment/atom
<compute_cluster_atom>`, :doc:`compute chunk/atom molecule
<compute_chunk_atom>`, or :doc:`fix evaporate <fix_evaporate>` can then
be applied for analysis or algorithmic purposes.

The *Nevery*, *Nrepeat*, and *Nfreq* values are used for an internally-
created :doc:`fix store/state <fix_store_state>` command that saves bond
history. If the fraction of times that is a bond is identified, out of
*Nrepeat* sampled timesteps, is greater than *fraction*, then this fix
assigns a bond. For example if *fraction* = 0.5, and a bond was
identified between atoms *i* and *j* for 60 of the 100 previously sampled
timesteps, then a bond is added between atom *i* and *j* on the current
timestep. The bonds identified by this fix are added as explicit bonds so
that the geometry of bond pairs is accessible to other commands. A bond
type of 1 is assigned to all bonds added by this fix. All existing bonds
that do not fulfill this criterion are removed.

.. note::

   This fix should generally not be used with models that require
   explicit bonds.  This fix adds and *deletes* bonds but does not
   update angles, dihedrals, or impropers.  This fix is meant to assist
   with on-the-fly analysis for atom styles that do not define explicit
   bonds, such as machine learning interatomic potentials or bond-order
   force fields. :doc:`Bond_style zero <bond_zero>` can be used to
   prevent these guessed bonds from affecting the dynamics of the
   simulation.

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
