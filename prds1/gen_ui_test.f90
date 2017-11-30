subroutine gen_ui_test( )
  !
  ! ui: ensemble of LHV
  ! u: distribution function over ui
  !
  use commvars, only : DP, nui, xstart, ndim, ui, pi
  !
  implicit none
  !
!  integer, intent(in) :: nui, xstart, ndim
!  real(DP), intent(out) :: ui(nui,xstart:ndim)
  !
  real(DP) :: mu(nui)
  real(DP) :: u, v, theta, phi
  integer :: i, j
  !
  allocate( ui(nui,xstart:ndim) )
  ui(:,:) = 0.d0
  !
  if ( xstart==1 ) then
    !
    ui(1,xstart:ndim) = (/ 1.d0, 1.d0, 0.d0 /)
    ui(2,xstart:ndim) = (/ 1.d0, 0.d0, 1.d0 /)
    ui(3,xstart:ndim) = (/ 1.d0,-1.d0, 0.d0 /)
    ui(4,xstart:ndim) = (/ 1.d0, 0.d0,-1.d0 /)
    !
    ui(:,xstart:ndim) = 0.25d0*ui(:,xstart:ndim)
    !
  elseif ( xstart==0 ) then
    !
    !ui(1,xstart:ndim) = (/ 1.d0,  1.d0, 0.d0, 0.d0 /)
    !ui(2,xstart:ndim) = (/ 1.d0,  0.d0, 1.d0, 0.d0 /)
    !ui(3,xstart:ndim) = (/ 1.d0,  0.d0, 0.d0, 1.d0 /)
    !ui(4,xstart:ndim) = (/ 1.d0,  -1.d0/sqrt(3.d0), -1.d0/sqrt(3.d0), -1.d0/sqrt(3.d0) /)
    ui(1,xstart:ndim) = (/ 1.d0,  -1.d0/sqrt(3.d0),  1.d0/sqrt(3.d0),  1.d0/sqrt(3.d0) /)
    ui(2,xstart:ndim) = (/ 1.d0,   1.d0/sqrt(3.d0), -1.d0/sqrt(3.d0),  1.d0/sqrt(3.d0) /)
    ui(3,xstart:ndim) = (/ 1.d0,   1.d0/sqrt(3.d0),  1.d0/sqrt(3.d0), -1.d0/sqrt(3.d0) /)
    ui(4,xstart:ndim) = (/ 1.d0,  -1.d0/sqrt(3.d0), -1.d0/sqrt(3.d0), -1.d0/sqrt(3.d0) /)
    !
    ui(:,:) = 0.25d0*ui(:,:)
    !
  endif
  !
  return
  !
end subroutine gen_ui_test
