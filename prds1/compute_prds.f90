subroutine compute_prds( )
  ! ui: ensemble of LHV
  ! 
  use commvars, only: DP, eps8, eps10, xstart, ndim, ndim1, nui, invrho, ui, lunidis, &
                      n_extrpts, extrpts_coor, combi
  !
  implicit none
  !
  logical :: lintersect, lnested
  integer :: i, j, k, l, n, p, r, np, np2, nc, nr, ii, jj
  integer, allocatable :: indx(:,:)
  real(DP) :: d1, d2
  real(DP) :: prds
  !real(DP) :: M(3,1:ndim), M4(3,0:ndim), normal(1:ndim1), vup(0:7,0:ndim), vdown(0:7,0:ndim)
  real(DP) :: M(3,1:ndim), M4(20,0:ndim), vu(0:ndim), vd(0:ndim)
  real(DP) :: normalU(0:ndim1), normalD(0:ndim1), normal(1:ndim1)
  real(DP) :: points_coor(xstart:ndim,xstart:ndim)
  real(DP), allocatable :: vup(:,:), vdown(:,:), vupA(:,:), vdownA(:,:), prds_global(:)
  !
  integer, external :: ckn
  !
  integer :: TID, NTHREADS, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
  !
  ! create array indx
  !
!  n = 0; nc = ckn(3,nui); allocate( indx(nc,3) )
!  do i = 1, nui
!    do j = i+1, nui
!      do k = j+1, nui
!        n = n + 1
!        indx(n,1) = i; indx(n,2) = j; indx(n,3) = k
!      enddo ! k loop
!    enddo ! j loop
!  enddo ! i loop
  !
  ! Calculate the extreme points
  !
!$OMP PARALLEL PRIVATE(NTHREADS, TID, i, j, k, p, np, np2, l, n, nr, ii, jj, r), &
!$OMP          PRIVATE(combi, lintersect, d1, d2, prds, points_coor, M, M4, normal), &
!$OMP          PRIVATE(normalU, normalD, vu, vd, vup, vdown, vupA, vdownA), &
!$OMP          SHARED( nui, ui, indx, prds_global, xstart ) 
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
!  do p = 1, nc

  do i = 1, nui
    do j = i+1, nui
      do k = j+1, nui
        !
!        M(1,1:ndim) = ui(indx(p,1),1:ndim); M4(1,0:ndim) = ui(indx(p,1),0:ndim)
!        M(2,1:ndim) = ui(indx(p,2),1:ndim); M4(2,0:ndim) = ui(indx(p,2),0:ndim)
!        M(3,1:ndim) = ui(indx(p,3),1:ndim); M4(3,0:ndim) = ui(indx(p,3),0:ndim)
        M(1,1:ndim) = ui(i,1:ndim); M4(1,0:ndim) = ui(i,0:ndim)
        M(2,1:ndim) = ui(j,1:ndim); M4(2,0:ndim) = ui(j,0:ndim)
        M(3,1:ndim) = ui(k,1:ndim); M4(3,0:ndim) = ui(k,0:ndim)
        !
        ! normal vector to the plane (M1M2M3) in the Bloch hyperplane
        call plane_normal_vect( M(1:3,1:ndim), normal(1:ndim1) )
        !
        ! compute 2 extreme points from upper part (UP) and lower part (DOWN)
        ! of the plane and count all points belonging to the plane  
        vu(:) = 0.d0; vd(:) = 0.d0
        np = 3 ! count points on the same plane
        do l = 1, nui
          !
!          if ( (l.eq.indx(p,1)) .or. (l.eq.indx(p,2)) .or. (l.eq.indx(p,3)) ) cycle
          if ( (l.eq.i) .or. (l.eq.j) .or. (l.eq.k) ) cycle
          d1 = dot_product( normal(1:ndim), ui(l,1:ndim) ) - normal(ndim1)
          if ( abs(d1).lt.eps10 ) then
!            write(*,*)"On the same plane", d1
            np = np + 1
            if ( np.gt.20 ) write(*,*)"np TOO LARGE: out-of-bound!"
            M4(np,0:ndim) = ui(l,0:ndim)
            !cycle
          else
            if ( d1.gt.0.d0) then
              vu(0) = vu(0) + ui(l,0) ! note: dmu(n) goes into ui
              vu(1:ndim) = vu(1:ndim) + ui(l,1:ndim) 
            else
              vd(0) = vd(0) + ui(l,0) 
              vd(1:ndim) = vd(1:ndim) + ui(l,1:ndim) 
            endif
          endif
          !
        enddo
        !
        ! generate new extreme points by adding M-point to vup and vdown
        !
        np2 = 2**np
        allocate( vup(0:np2-1,0:ndim), vdown(0:np2-1,0:ndim) )
        vup(:,:) = 0.d0; vdown(:,:) = 0.d0
        vup(0,0:ndim) = vu(0:ndim); vdown(0,0:ndim) = vd(0:ndim) 
        !allocate( vup(0:ndim,0:2**np-1), vdown(0:ndim,0:2**np-1) )
        !vup(:,:) = 0.d0; vdown(:,:) = 0.d0
        !vup(0:ndim,0) = vu(0:ndim); vdown(0:ndim,0) = vd(0:ndim) 
        !
        n = 0
        do r = 1, np
           !
           nr = ckn( r, np ); allocate( combi(nr,r))
!write(*,*)r, np, nr
           call combinations( r, np, nr, combi )
           !
           do ii = 1, nr
             n = n + 1
               vup(n,0:ndim) =   vup(0,0:ndim) 
             vdown(n,0:ndim) = vdown(0,0:ndim) 
             !  vup(0:ndim,n) =   vup(0:ndim,0) 
             !vdown(0:ndim,n) = vdown(0:ndim,0) 
             do jj = 1, r
                 vup(n,0:ndim) =   vup(n,0:ndim) + M4(combi(ii,jj),0:ndim)
               vdown(n,0:ndim) = vdown(n,0:ndim) + M4(combi(ii,jj),0:ndim)
              !   vup(0:ndim,n) =   vup(0:ndim,n) + M4(combi(ii,jj),0:ndim)
              !vdown(0:ndim,n) = vdown(0:ndim,n) + M4(combi(ii,jj),0:ndim)
             enddo
           enddo
           !
           deallocate( combi )
           !
        enddo 
        !
        ! bring these points to Alice's space
        !
        allocate( vupA(0:ndim,0:np2-1), vdownA(0:ndim,0:np2-1) )
        vupA(:,:) = 0.d0; vdownA(:,:) = 0.d0
        call dgemm( 'N', 'T', 4, np2, 4, 1.d0, invrho(0:ndim,0:ndim), 4, &
                     vup(0:np2-1,0:ndim), np2, 0.d0, vupA(0:ndim,0:np2-1),   4 )
        !call dgemm( 'N', 'N', 4, 2**np, 4, 1.d0, invrho(0:ndim,0:ndim), 4, &
        !             vup(0:ndim,0:2**np-1),   4, 0.d0, vupA(0:ndim,0:2**np-1),   4 )
        !vupA(0:ndim,0:2**np-1) = matmul( invrho(0:ndim,0:ndim), transpose(vup(0:2**np-1,0:ndim)) )
        call dgemm( 'N', 'T', 4, np2, 4, 1.d0, invrho(0:ndim,0:ndim), 4, &
                     vdown(0:np2-1,0:ndim), np2, 0.d0, vdownA(0:ndim,0:np2-1), 4 )
        !call dgemm( 'N', 'N', 4, 2**np, 4, 1.d0, invrho(0:ndim,0:ndim), 4, &
        !             vdown(0:ndim,0:2**np-1), 4, 0.d0, vdownA(0:ndim,0:2**np-1), 4 )
        !vdownA(0:ndim,0:2**np-1) = matmul( invrho(0:ndim,0:ndim), transpose(vdown(0:2**np-1,0:ndim)) )
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
        ! check if "DOWN" facet intersects with the Bloch hyperplane 
        !
        lintersect = .false.
        do l = 0, np2-1
          d1 = vdownA(0,l) - 1.d0; if ( abs(d1).gt.eps8 ) exit
        enddo
        do l = 0, np2-1
          if ( abs(vdownA(0,l)-1.d0).gt.eps8 ) d2 = vdownA(0,l) - 1.d0
          if ( d1*d2.lt.0.d0 ) then
            lintersect = .true.;  exit
          endif
        enddo
        !
        if ( lintersect ) then
          !
          ! equation for the "UP" hyperplane
          do l = 0, np2-4
            !
            points_coor(xstart,  xstart:ndim) = vdownA(xstart:ndim,l) 
            points_coor(xstart+1,xstart:ndim) = vdownA(xstart:ndim,l+1) 
            points_coor(xstart+2,xstart:ndim) = vdownA(xstart:ndim,l+2) 
            points_coor(xstart+3,xstart:ndim) = vdownA(xstart:ndim,l+3) 
            !
            ! compute normal vector of the hyperplane and the right hand side
            call hyp_normal_vect( points_coor(xstart:ndim,xstart:ndim), normalD(xstart:ndim1) )
            !
            ! if the normal vector is NOT zero, i.e. these points are  enough to determine
            ! a hyperplane, exit the do loop
            d2 = dot_product( normalD(xstart:ndim), normalD(xstart:ndim) )
            if ( abs(d2).gt.eps10 ) exit
            !
          enddo
          !
          ! distance from (1,0,0,0) to the 3D plane
          d1 = abs(normalD(xstart)-normalD(ndim1))/sqrt(dot_product( normalD(1:ndim),normalD(1:ndim)))
          if ( d1.lt.prds ) prds = d1
          !
        endif
        !
        deallocate( vup, vdown, vupA, vdownA )
        !
!  enddo ! p loop
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
!  deallocate( indx )
  !
  write(*,'(A,I5,f18.14)') "PRINCIPAL RADIUS: ", nui, minval(prds_global(0:))
  !
end subroutine compute_prds
