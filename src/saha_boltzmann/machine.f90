!> Machine-dependent constants.
MODULE machine
  IMPLICIT NONE
  INTEGER, PARAMETER :: sp = selected_real_kind(precision(1.0)) !< single precision
  INTEGER, PARAMETER :: dp = selected_real_kind(precision(1.0d0)) !< double precision
END MODULE machine
