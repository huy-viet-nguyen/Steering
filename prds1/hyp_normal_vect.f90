subroutine plane_normal_vect( points_coor, normal )
  !
  use commvars,  only : DP, ndim, ndim1
  !
  implicit none
  !
  real(DP), intent(in) :: points_coor(1:ndim,1:ndim)
  real(DP), intent(out) :: normal(1:ndim1)
  ! last element is used for the offset
  !
  real(DP) :: v(2,1:ndim)
  integer :: i, j, k
  !
  ! two vectors in the plane
  v(1,1:ndim) = points_coor(1,1:ndim)-points_coor(ndim,1:ndim)
  v(2,1:ndim) = points_coor(2,1:ndim)-points_coor(ndim,1:ndim)
  !
  ! cross product
  normal(1) = v(1,2)*v(2,3) - v(1,3)*v(2,2)
  normal(2) = v(1,3)*v(2,1) - v(1,1)*v(2,3)
  normal(3) = v(1,1)*v(2,2) - v(1,2)*v(2,1)
  normal(1:3) = normal(1:3)/sqrt(dot_product(normal(1:3),normal(1:3)))
  !
  normal(ndim1) = dot_product(normal(1:ndim),points_coor(1,1:ndim))
  !
end subroutine plane_normal_vect
!
subroutine hyp_normal_vect( points_coor, normal )
  !
  use commvars,  only : DP, ndim, ndim1, xstart
  !
  implicit none
  !
  real(DP), intent(in) :: points_coor(xstart:ndim,xstart:ndim)
  real(DP), intent(out) :: normal(xstart:ndim1)
  ! last element is used for the offset
  real(DP), external :: det3
  !
  real(DP) :: v1(xstart:ndim), v2(xstart:ndim), v3(xstart:ndim)
  real(DP) :: v(3,xstart:ndim)
  integer :: i, j, k
  !
  !
  v(1,0:ndim) = points_coor(1,0:ndim)-points_coor(0,0:ndim)
  v(2,0:ndim) = points_coor(2,0:ndim)-points_coor(0,0:ndim)
  v(3,0:ndim) = points_coor(3,0:ndim)-points_coor(0,0:ndim)
  !
  normal(0) =  det3( v(:,1), v(:,2), v(:,3) )
  normal(1) = -det3( v(:,0), v(:,2), v(:,3) )
  normal(2) =  det3( v(:,0), v(:,1), v(:,3) )
  normal(3) = -det3( v(:,0), v(:,1), v(:,2) )
  normal(0:3) = normal(0:3)/sqrt(dot_product(normal(0:3),normal(0:3)))
  !
  normal(ndim1) = dot_product(normal(0:ndim),points_coor(1,0:ndim))
  !
  !
end subroutine hyp_normal_vect
!
!
!
function det3 ( A1, A2, A3 )
  ! A1, A2, A3: 3 column vectors that forms a 3x3-matrix A
  ! det: determinant of A
  !
  use commvars,  only: DP
  !
  implicit none
  !
  real(DP), intent(in) :: A1(3), A2(3), A3(3)
  real(DP) :: det3
  !
  det3 = A1(1)*A2(2)*A3(3) + &
         A1(2)*A2(3)*A3(1) + &
         A1(3)*A2(1)*A3(2) - &
         A3(1)*A2(2)*A1(3) - &
         A3(2)*A2(3)*A1(1) - &
         A3(3)*A2(1)*A1(2)
  !
end function det3
