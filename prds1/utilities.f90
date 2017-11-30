subroutine combinations( r, n, crn, combi ) 
  !
  implicit none
  !
  integer, intent(in) :: r, n, crn
  integer, intent(out) :: combi(crn,r)
  integer :: np, m, ii
  integer, allocatable :: comb(:)     
  !
  m = 1
  do ii = 0, r-1
    m = m*(n-ii)/(ii+1)
  enddo
  if ( m.ne.crn ) then
    write(*,*) "crn is not equal to C^r_n"
    stop
  endif
  !
  allocate(comb(r))     
  !
  np = 0
  call iterate(1,n-r+1,1) 
  !
  contains 
    !
    recursive subroutine iterate(s,e,j) 
      !
      integer, intent(in) :: s, e, j       
      integer :: i
      !
      do i = s, e          
        comb(j) = i          
        if ( j.lt.r ) call iterate( i+1, e+1, j+1 )      
        if ( j.eq.r ) then
          !WRITE(*,*)comb      
          np = np + 1
          combi(np,:) = comb(:)
        endif
      end do 
    end subroutine iterate 
end subroutine combinations
!
!
!
integer function ckn( k_, n_ )
  !
  implicit none
  ! 
  integer, intent(in) :: k_, n_
  !
  integer :: i_, m_
  !
  m_ = 1
  do i_ = 0, k_-1
    m_ = m_*(n_-i_)/(i_+1)
  enddo
  !
  ckn = m_
  !
end function ckn
