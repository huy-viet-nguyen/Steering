subroutine gen_ui_uniform( )
  !
  ! ui: ensemble of LHV
  ! u: distribution function over ui
  !
  use commvars, only : DP, eps8, nui, xstart, ndim, ui, pi
  !
  implicit none
  !
  real(DP) :: mu(nui)
  real(DP) :: u, v, theta, phi
  integer :: i, j
  !
  allocate( ui(nui,xstart:ndim) )
  ui(:,:) = 0.d0
  !
  ! generate points of uniform distribution on a 3D unit sphere
  !
  !if ( mod(nui,2).ne.0 ) then
  !  write(*,*) "nui must be even"
  !  stop
  !endif
  !
  ui(:,0) = 1.d0
  do i = 1, nui
    !
    u = randy(); v = randy()
    phi = 2.d0*pi*u
    theta = acos(2.d0*v-1.d0)
    if ( (2.d0*theta).gt.pi) theta = pi - theta
    !
    ui(i,1) = cos(theta) 
    ui(i,2) = sin(theta)*cos(phi)
    ui(i,3) = sin(theta)*sin(phi)
    !
    !ui(i+nui/2,1) = -cos(theta)
    !ui(i+nui/2,2) = -sin(theta)*cos(phi) 
    !ui(i+nui/2,3) = -sin(theta)*sin(phi)
    !
  enddo
  !
  ! check to make sure that all the points are on Bob's positive cone 
  !
  do i = 1, nui
    u = ui(i,0)**2 - dot_product( ui(i,1:ndim),ui(i,1:ndim) )
    if ( abs(u).gt.eps8 ) then
      write(*,*)"Not on Bob's positive cone", i
      stop
    endif
  enddo
  !
  mu(:) = 1.d0/dble(2*nui)
  do i = 1, nui
    ui(i,:) = mu(i)*ui(i,:)
  enddo
  !
contains
    !
    !------------------------------------------------------------------------
    REAL(DP) FUNCTION randy ( irand )
      !------------------------------------------------------------------------
      !
      ! x=randy(n): reseed with initial seed idum=n ( 0 <= n <= ic, see below)
      !             if randy is not explicitly initialized, it will be
      !             initialized with seed idum=0 the first time it is called
      ! x=randy() : generate uniform real(DP) numbers x in [0,1]
      !
      implicit none
      !
      integer, parameter :: dp = selected_real_kind(14,200)
      !REAL(DP) :: randy
      INTEGER, optional    :: irand
      !
      INTEGER , PARAMETER  :: m    = 714025, &
                              ia   = 1366, &
                              ic   = 150889, &
                              ntab = 97
      REAL(DP), PARAMETER  :: rm = 1.0_DP / m
      INTEGER              :: j
      INTEGER, SAVE        :: ir(ntab), iy, idum=0
      LOGICAL, SAVE        :: first=.true.
      !
      IF ( present(irand) ) THEN
         idum = MIN( ABS(irand), ic)
         first=.true.
      ENDIF

      IF ( first ) THEN
         !
         first = .false.
         idum = MOD( ic - idum, m )
         !
         DO j=1,ntab
            idum=mod(ia*idum+ic,m)
            ir(j)=idum
         END DO
         idum=mod(ia*idum+ic,m)
         iy=idum
      END IF
      j=1+(ntab*iy)/m
      !IF( j > ntab .OR. j <  1 ) call errore('randy','j out of range',ABS(j)+1)
      IF( j > ntab .OR. j <  1 ) write(*,*) "randy: j out of range"
      iy=ir(j)
      randy=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
      !
      RETURN
      !
    END FUNCTION randy
    !
end subroutine gen_ui_uniform
