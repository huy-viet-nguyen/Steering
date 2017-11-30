subroutine setup_rho( )
  !
  use commvars,   only : DP, rho, invrho
  !
  implicit none
  !
  real(DP) :: p
  !
  p = 0.49d0
  !
  rho(:,:) = 0.d0
  invrho(:,:) = 0.d0
  !
  rho(0,0) = 0.5d0
  rho(1,1) = 0.5d0*p
  rho(2,2) =-0.5d0*p
  rho(3,3) = 0.5d0*p
  !
  invrho(0,0) = 1.d0/rho(0,0)
  invrho(1,1) = 1.d0/rho(1,1)
  invrho(2,2) = 1.d0/rho(2,2)
  invrho(3,3) = 1.d0/rho(3,3)
  !
end subroutine setup_rho
