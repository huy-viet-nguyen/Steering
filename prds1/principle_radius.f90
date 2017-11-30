program extrpts_gen
  !
  use commvars, only : DP, xstart, nui, ndim, ndim1, combi, &
                       lunidis, lnew, lblas
  !
  implicit none
  !
  integer :: i
  !
  ! input readin
  !
  read(5,*) nui
  read(5,*) lunidis
  read(5,*) lnew
  read(5,*) lblas
  !
  if ( lblas.and.(nui.lt.3) ) then
    write(*,*) "WARNING: nui too small for BLAS routine"
    stop
  endif
  !
  ! check combinations routine
  !
  !allocate ( combi(10,3))
  !call combinations( 3, 5, 10, combi)
  !do i = 1, 10
  !  write(*,*)combi(i,:)
  !enddo
  !stop
  !
  ! Calculation in 4D
  !
  xstart = 0 
  !
  call setup_rho( )
  !
!  if ( nui.eq.4 ) then
!    call gen_ui_test( )
!    call box_u_bf ( )
!  endif
  !
  if ( lunidis ) then
    call gen_ui_uniform( )
  else
    call gen_ui_lebedev( )
  endif
  !call compute_extreme_points( )
  !call compute_extreme_points_new( )
  call compute_prds_half( )
  !call box_u ( )
  !
  !
end program extrpts_gen
