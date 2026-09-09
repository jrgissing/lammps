.. index:: compute bond

compute guess/bonds command
====================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID guess/bonds mode values

* ID, group-ID are documented in :doc:`compute <compute>` command
* guess/bonds = style name of this compute command
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

   compute 1 all guess/bonds distance 1.15 c1 0.70 hc 0.37 oc 0.6 n 0.65
   compute 1 all guess/bonds pauling 0.3 C C 1.54 0.353 C H 1.09 0.399 C O 1.43 0.279 C N 1.47 0.298 O N 1.44 0.26

Description
"""""""""""

Guess whether or not a covalent bond exists between each atom pair. This
compute outputs an array of the calculated bond connectivity. It does
not actually assign bonds to the bond list, `fix assign/bonds
<fix_assign_bonds>` can be used to do add the explicit bonds identified
by this fix.

In *distance* mode, if a bond satisfies the distance cutoff criterion on
a given timestep, it is a bond is guessed to exist for that atom pair. In
*distance* mode, the distance cutoff is defined by the following
equation:

.. math::

   cutoff = prefactor*(r_1 + r_2)

where :math:`r_1` is the atomic radius of the first atom and :math:`r_2`
is the atomic radius of the second atom, and a bond is listed in the
compute output between these two atoms if their interatomic distance is
less than the cutoff. The prefactor must be provided directly after the
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
:math:`b` (see example above). An bond between two atoms is listed in the
compute if their calculated bond order is less than *bond_order_cutoff*.


Output info
"""""""""""

This compute calculates a per-atom array with :math:`N` columns, where
:math:`N` is the maximum number of bonds per atom plus one. The first
column is the number of bonds guessed for that atom, while columns 2-N
are the atom IDs of atoms guessed to be bonded to that atom. These values
can be accessed by indices 1 through :math:`N` by any command that uses
per-atom values from a compute as input.  See the :doc:`Howto output
<Howto_output>` page for an overview of LAMMPS output options.

Restrictions
""""""""""""
none

Related commands
""""""""""""""""

:doc:`fix assign/bonds <fix_assign_bonds>`

Default
"""""""

none
