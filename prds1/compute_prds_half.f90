subroutine compute_prds_half( )
  ! ui: ensemble of LHV
  ! 
  use commvars, only: DP, eps8, eps10, xstart, ndim, ndim1, nui, invrho, ui, lunidis, &
                      n_extrpts, extrpts_coor, combi
  !
  implicit none
  !
  real(DP) :: d1, d2
  real(DP) :: prds
  real(DP) :: M(3,1:ndim), vu(0:ndim)
  real(DP) :: normalU(0:ndim1), normal(1:ndim1)
  real(DP) :: points_coor(xstart:ndim,xstart:ndim)
  real(DP), allocatable :: vup(:,:), vupA(:,:), prds_global(:)
  !
  integer :: i, j, k, l, n, p, r, np, np2, nc, nr, ii, jj
  integer :: indx(90)
  integer, external :: ckn
  integer :: TID, NTHREADS, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
  !
  logical :: lintersect
  !
  ! Calculate the extreme points
  !
!$OMP PARALLEL PRIVATE(NTHREADS, TID, i, j, k, p, np, np2, l, n, nr, ii, jj, r), &
!$OMP          PRIVATE(combi, lintersect, d1, d2, prds, points_coor, M, normal), &
!$OMP          PRIVATE(normalU, vu, vup, vupA, indx ), &
!$OMP          SHARED( nui, ui, prds_global, xstart ) 
  !
  TID = OMP_GET_THREAD_NUM()
  WRITE(*,*) "THREAD ID = ", TID
  IF (TID .EQ. 0) THEN
    NTHREADS = OMP_GET_NUM_THREADS()
    ALLOCATE( prds_global(0:NTHREADS) )
    prds_global(0:NTHREADS) = 10000000.d0
  ENDIF
  prds = 10000000.d0
  !
!$OMP DO 
  !
  do i = 1, nui
    do j = i+1, nui
      do k = j+1, nui
        !
        M(1,1:ndim) = ui(i,1:ndim); indx(1) = i
        M(2,1:ndim) = ui(j,1:ndim); indx(2) = j
        M(3,1:ndim) = ui(k,1:ndim); indx(3) = k
        !
        ! normal vector to the plane (M1M2M3) in the Bloch hyperplane
        call plane_normal_vect( M(1:3,1:ndim), normal(1:ndim1) )
        !
        ! compute the extreme point from upper part (UP) of the sphere
        ! and count all points belonging to the plane  
        vu(:) = 0.d0; np = 3 ! count points on the same plane
        do l = 1, nui
          !
          if ( (l.eq.i) .or. (l.eq.j) .or. (l.eq.k) ) cycle
          d1 = dot_product( normal(1:ndim), ui(l,1:ndim) ) - normal(ndim1)
          if ( abs(d1).lt.eps10 ) then
!            write(*,*)"On the same plane", d1
            np = np + 1; indx(np) = l
            if ( np.gt.90 ) write(*,*)"np TOO LARGE: indx array out-of-bound!"
          else
            if ( d1.gt.0.d0) then
              vu(0) = vu(0) + ui(l,0) ! note: dmu(n) goes into ui
              vu(1:ndim) = vu(1:ndim) + ui(l,1:ndim) 
            endif
          endif
          !
        enddo
        !
        ! generate new extreme points by adding M-point to vup and vdown
        !
        np2 = 2**np
        allocate( vup(0:np2-1,0:ndim) ) 
        vup(:,:) = 0.d0; vup(0,0:ndim) = vu(0:ndim)
        !
        n = 0
        do r = 1, np
           !
           nr = ckn( r, np ); allocate( combi(nr,r))
           call combinations( r, np, nr, combi )
           !
           do ii = 1, nr
             n = n + 1
               vup(n,0:ndim) =   vup(0,0:ndim) 
             do jj = 1, r
                 !vup(n,0:ndim) =   vup(n,0:ndim) + M4(combi(ii,jj),0:ndim)
                 vup(n,0:ndim) =   vup(n,0:ndim) + ui(indx(combi(ii,jj)),0:ndim)
             enddo
           enddo
           !
           deallocate( combi )
           !
        enddo 
        !
        ! bring these points to Alice's space
        !
        allocate( vupA(0:ndim,0:np2-1) )
        vupA(:,:) = 0.d0
        call dgemm( 'N', 'T', 4, np2, 4, 1.d0, invrho(0:ndim,0:ndim), 4, &
                     vup(0:np2-1,0:ndim), np2, 0.d0, vupA(0:ndim,0:np2-1),   4 )
        !
        ! check if "UP" facet intersects with the Bloch hyperplane 
        !
        lintersect = .false.
        do l = 0, np2-1
          d1 = vupA(0,l) - 1.d0; if ( abs(d1).gt.eps8 ) exit
        enddo
        do l = 0, np2-1
          if ( abs(vupA(0,l)-1.d0).gt.eps8 ) d2 = vupA(0,l) - 1.d0
          if ( d1*d2.lt.0.d0 ) then
            lintersect = .true.;  exit
          endif
        enddo
        !
        ! if YES, compute the principle radius
        !
        if ( lintersect ) then
          !
          ! equation for the "UP" hyperplane. the do loop is to make sure that 
          ! we find 4 points that are not on a 3D plane
          do l = 0, np2-4
            !
            points_coor(xstart,  xstart:ndim) = vupA(xstart:ndim,l) 
            points_coor(xstart+1,xstart:ndim) = vupA(xstart:ndim,l+1) 
            points_coor(xstart+2,xstart:ndim) = vupA(xstart:ndim,l+2) 
            points_coor(xstart+3,xstart:ndim) = vupA(xstart:ndim,l+3) 
            !
            ! compute normal vector of the hyperplane and the right hand side
            call hyp_normal_vect( points_coor(xstart:ndim,xstart:ndim), normalU(xstart:ndim1) )
            !
            ! if the normal vector is NOT zero, i.e. these points are  enough to determine
            ! a hyperplane, exit the do loop
            d2 = dot_product( normalU(xstart:ndim), normalU(xstart:ndim) )
            if ( abs(d2).gt.eps10 ) exit
            !
          enddo
          !
          ! distance from (1,0,0,0) to the 3D plane
          d1 = abs(normalU(xstart)-normalU(ndim1))/sqrt(dot_product( normalU(1:ndim),normalU(1:ndim)))
          if ( d1.lt.prds ) prds = d1
          !
        endif
        !
        deallocate( vup, vupA )
        !
      enddo ! k loop
    enddo ! j loop
  enddo ! i loop
  !
!$OMP END DO
  ! 
  prds_global(TID) = prds
write(*,*) "prds =", prds
  !
!$OMP END PARALLEL
  !
  write(*,'(A,I5,f18.14)') "PRINCIPAL RADIUS (HALF): ", nui, minval(prds_global(0:))
  !
end subroutine compute_prds_half
