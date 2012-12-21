module precision_mod
implicit none
  INTEGER, PARAMETER :: sp = selected_real_kind(precision(1.0)) !< single precision
  INTEGER, PARAMETER :: dp = selected_real_kind(precision(1.0d0)) !< double precision
end module precision_mod
