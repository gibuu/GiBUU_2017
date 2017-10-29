program test
  use inputGeneral, only: readInputGeneral
  use pi0eta_photoproduction
  implicit none
  real :: W

  call readInputGeneral()

  open(10, file="pi0eta.dat")
  write(10,'(2A15)') "# W[GeV] "  , "sigma [mu b]" 
  W=0.
  do
     W=W+0.01
     write(10,'(2E15.5)') W, sigma_gamma_p_to_p_pi0_eta(W)
     if (W.gt.4.4) exit
  end do
  close(10)
end program test
