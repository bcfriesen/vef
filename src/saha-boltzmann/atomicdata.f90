!> Holds all the relevant atomic data and physical constants. All units are CGS
!! unless specified explicitly.
MODULE atomicdata
  USE precision_mod
  IMPLICIT NONE

  INTEGER, PARAMETER :: nlvlmax = 100 !< maximum number of levels for any given atom.

  INTEGER :: nlvl(92, 0:92) !< actual number of levels to consider for ion \f$ j \f$ of
                            !! atom \f$ i \f$. nlvl\f$ (i,j) = 0 \f$ if \f$ j > i \f$. \f$ n \f$ is equal to the
                            !! index \f$ k \f$ in the formalism from class. Just for future reference,
                            !! \f$ i = 1 \rightarrow \f$ hydrogen, \f$ j = 0 \rightarrow \f$ neutral,
                            !! \f$ k = 1 \rightarrow \f$ ground state.
                            !! hence the non-standard indices for \f$ j \f$ in a whole bunch of these
                            !! arrays.
  REAL (KIND=dp) :: chilvl(92, 0:92, nlvlmax) !< level energies (ground = 0) (units = eV)
  INTEGER :: g_i(92, 0:92, nlvlmax) !< statistical weight for each level
  REAL (KIND=dp) :: chiion(92, 0:92) !< chiion\f$(i,j)\f$: ionization potential of atom \f$i\f$ in
                                     !! ionization state \f$j\f$ (units: eV)
  INTEGER, PARAMETER :: g_e = 2 !< electron statistical weight
END MODULE atomicdata
