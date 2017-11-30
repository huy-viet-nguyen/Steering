subroutine gen_ui_lebedev( )
  !
  ! ui: ensemble of LHV
  ! u: distribution function over ui
  !
  use commvars, only : DP, nui, eps8, eps10, xstart, ndim, ndim1, ui, pi
  !
  implicit none
  !
  real(DP), allocatable :: mu(:), x(:), y(:), z(:), w(:), plane(:,:)
  real(DP) :: M(1:3,1:ndim), M4(1:20,0:ndim), normal(1:ndim1)
  real(DP) :: d1
  logical :: lrepeat
  integer :: i, j, k, l, n, np, npl, nz0
  integer, external :: ckn
  !
  ! generate Lebedev abscissa and weights
  !
  write(*,*) "UNIFORM DISTRIBUTION MODELED BY LEBEDEV ABSCISSA"
  !
  allocate( mu(nui), x(nui), y(nui), z(nui), w(nui) )
  call ld_by_order ( nui, x, y, z, w )
  mu(:) = 1.d0/dble(nui)
  !
  allocate( ui(nui,xstart:ndim) )
  ui(:,:) = 0.d0
  !
  do i = 1, nui
    ui(i,0) = mu(i); ui(i,1) = z(i)*mu(i)
    ui(i,2) = x(i)*mu(i);  ui(i,3) = y(i)*mu(i)
  enddo
  !
  deallocate( mu, x, y, z, w )
  !
if ( .false. ) then 
  ! check coplanar points
  !
  allocate( plane(ckn(3,nui),1:4) )
  npl = 0
  !
  do i = 1, nui
    do j = i+1, nui
      do k = j+1, nui
        !
        M(1,1:ndim) = ui(i,1:ndim); M4(1,0:ndim) = ui(i,0:ndim)
        M(2,1:ndim) = ui(j,1:ndim); M4(2,0:ndim) = ui(j,0:ndim)
        M(3,1:ndim) = ui(k,1:ndim); M4(3,0:ndim) = ui(k,0:ndim)
        !
        ! normal vector to the plane (M1M2M3) in the Bloch hyperplane
        call plane_normal_vect( M(1:3,1:ndim), normal(1:ndim1) )
        !
        ! count all points belonging to the same plane  
        np = 3 ! count points on the same plane
        do l = 1, nui
          !
          if ( (l.eq.i) .or. (l.eq.j) .or. (l.eq.k) ) cycle
          !
          d1 = dot_product( normal(1:ndim), ui(l,1:ndim) ) - normal(ndim1)
          if ( abs(d1).lt.eps10 ) then
!            write(*,*)"On the same plane", d1
            np = np + 1
            if ( np.gt.20 ) write(*,*)"np TOO LARGE: out-of-bound!"
            M4(np,0:ndim) = ui(l,0:ndim)
          endif
          !
        enddo
        !
        if ( np.gt.3 ) then
          !
          lrepeat = .false.
          do l = 1, npl
            d1 = dot_product( normal(1:3), plane(l,1:3) )
            if ( (abs(1.d0-d1).lt.eps8).and.(abs(normal(4)-plane(l,4)).lt.eps8 ) .or. &
                 (abs(1.d0+d1).lt.eps8).and.(abs(normal(4)+plane(l,4)).lt.eps8 ) ) then
              lrepeat = .true.
              exit
            endif 
          enddo
          !
          if ( .not.lrepeat ) then
            npl = npl + 1
            plane(npl,1:4) = normal(1:4)
            !
            write(*,*) "The number of points on this plane: ", np
            write(*,'(A,4(f16.12))')"normal", normal(1:ndim1)
            do n = 1, np
              write(*,'(4(f16.12))') M4(n,0:ndim)
            enddo
          endif
          !
        endif
        !
      enddo
    enddo
  enddo
  !
  deallocate( plane )
endif
  !
  return
  !
end subroutine gen_ui_lebedev
!
!
subroutine gen_ui_lebedev_invsym( )
  !
  ! ui: ensemble of LHV
  ! u: distribution function over ui
  !
  use commvars, only : DP, nui, eps8, eps10, xstart, ndim, ndim1, ui, pi
  !
  implicit none
  !
  real(DP), allocatable :: mu(:), x(:), y(:), z(:), w(:), plane(:,:)
  real(DP) :: M(1:3,1:ndim), M4(1:20,0:ndim), normal(1:ndim1)
  real(DP) :: d1
  logical :: lrepeat
  integer :: i, j, k, l, n, np, npl, nz0
  integer, external :: ckn
  !
  ! generate Lebedev abscissa and weights
  !
  write(*,*) "UNIFORM DISTRIBUTION MODELED BY LEBEDEV ABSCISSA"
  !
  allocate( mu(nui), x(nui), y(nui), z(nui), w(nui) )
  call ld_by_order ( nui, x, y, z, w )
  mu(:) = 1.d0/dble(nui)
  !
  ! np: number of points on z=0 plane
  !
  nz0 = count(abs(z).lt.eps10)
  !write(*,*)"nz0 = ", nz0
  allocate( ui((nui+nz0)/2,xstart:ndim) )
  ui(:,:) = 0.d0
  !
  np = 0
  do i = 1, nui
  !write(*,*)x(i), y(i), z(i)
    if ( z(i).gt.(-eps10) ) then
      np = np + 1
      ui(np,0) = mu(i); ui(np,1) = z(i)
      ui(np,2) = x(i);  ui(np,3) = y(i)
    endif
  enddo
  if ( np.ne.(nui+nz0)/2 ) then
    write(*,*)" Lebedev: np /= nui", np, (nui+nz0)/2
    stop
  else
    nui = (nui+nz0)/2
  endif
  !
  deallocate( mu, x, y, z, w )
  !
if ( .false. ) then 
  ! check coplanar points
  !
  allocate( plane(ckn(3,nui),1:4) )
  npl = 0
  !
  do i = 1, nui
    do j = i+1, nui
      do k = j+1, nui
        !
        M(1,1:ndim) = ui(i,1:ndim); M4(1,0:ndim) = ui(i,0:ndim)
        M(2,1:ndim) = ui(j,1:ndim); M4(2,0:ndim) = ui(j,0:ndim)
        M(3,1:ndim) = ui(k,1:ndim); M4(3,0:ndim) = ui(k,0:ndim)
        !
        ! normal vector to the plane (M1M2M3) in the Bloch hyperplane
        call plane_normal_vect( M(1:3,1:ndim), normal(1:ndim1) )
        !
        ! count all points belonging to the same plane  
        np = 3 ! count points on the same plane
        do l = 1, nui
          !
          if ( (l.eq.i) .or. (l.eq.j) .or. (l.eq.k) ) cycle
          !
          d1 = dot_product( normal(1:ndim), ui(l,1:ndim) ) - normal(ndim1)
          if ( abs(d1).lt.eps10 ) then
!            write(*,*)"On the same plane", d1
            np = np + 1
            if ( np.gt.20 ) write(*,*)"np TOO LARGE: out-of-bound!"
            M4(np,0:ndim) = ui(l,0:ndim)
          endif
          !
        enddo
        !
        if ( np.gt.3 ) then
          !
          lrepeat = .false.
          do l = 1, npl
            d1 = dot_product( normal(1:3), plane(l,1:3) )
            if ( (abs(1.d0-d1).lt.eps8).and.(abs(normal(4)-plane(l,4)).lt.eps8 ) .or. &
                 (abs(1.d0+d1).lt.eps8).and.(abs(normal(4)+plane(l,4)).lt.eps8 ) ) then
              lrepeat = .true.
              exit
            endif 
          enddo
          !
          if ( .not.lrepeat ) then
            npl = npl + 1
            plane(npl,1:4) = normal(1:4)
            !
            write(*,*) "The number of points on this plane: ", np
            write(*,'(A,4(f16.12))')"normal", normal(1:ndim1)
            do n = 1, np
              write(*,'(4(f16.12))') M4(n,0:ndim)
            enddo
          endif
          !
        endif
        !
      enddo
    enddo
  enddo
  !
  deallocate( plane )
endif
  !
  return
  !
end subroutine gen_ui_lebedev_invsym
