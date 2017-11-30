module commvars
!
implicit none
!
save
!
integer, parameter :: DP = selected_real_kind(14,200)
! integer, parameter :: ndim = 3, ndim1 = ndim+1, xstart = 0
integer, parameter :: ndim = 3, ndim1 = ndim+1
integer :: xstart
!
integer, allocatable :: combi(:,:)
!real(DP), parameter :: alp3(8,3) = transpose( reshape(        &
real(DP), parameter :: alp3(8,3) = transpose( reshape( (/ &
                          0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, &
                          0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, &
                          1.d0, 1.d0, 0.d0, 1.d0, 0.d0, 1.d0, &
                          0.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0 /), (/ 3, 8 /) ) )
real(DP), parameter :: alp4(16,4) = transpose( reshape( (/      &
                1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, &
                0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, &
                1.d0, 1.d0, 0.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0, &
                1.d0, 0.d0, 0.d0, 1.d0, 0.d0, 1.d0, 1.d0, 0.d0, &
                0.d0, 1.d0, 0.d0, 1.d0, 0.d0, 0.d0, 1.d0, 1.d0, &
                1.d0, 1.d0, 1.d0, 0.d0, 1.d0, 1.d0, 0.d0, 1.d0, &
                1.d0, 0.d0, 1.d0, 1.d0, 0.d0, 1.d0, 1.d0, 1.d0, &
                1.d0, 1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0 /), (/ 4, 16 /) ) )
!real(DP), parameter :: alp5(32,5) = 1.d0
!
real(DP), parameter :: eps(14) = (/ 1.d-1, 1.d-2, 1.d-3, 1.d-4, 1.d-5, 1.d-6, 1.d-7, 1.d-8, &
                                    1.d-9, 1.d-10, 1.d-11, 1.d-12, 1.d-13, 1.d-14 /) 
real(DP), parameter :: eps6 = 1.d-6, eps8 = 1.d-8, eps10 = 1.d-10
real(DP), parameter :: pi = acos(-1.d0)
!
integer :: nui = 5
logical :: lunidis=.true. ! uniform distribution LHS
logical :: lnew = .true. ! new scheme utilizing the symmetry
logical :: lblas = .true. ! exploit BLAS routine
!
real(DP) :: rho(0:ndim,0:ndim), invrho(0:ndim,0:ndim)
real(DP), allocatable :: ui(:,:)
!
integer :: n_extrpts
real(DP), allocatable :: extrpts_coor(:,:)
!
integer :: n_vertices, n_facets
real(DP), allocatable :: vertices_coor(:,:)
real(DP), allocatable :: facets_coor(:,:,:)
!
complex(DP), parameter :: one  = cmplx(1.d0,0.d0), &
                          zero = cmplx(0.d0,0.d0)
!
public :: DP, ndim, ndim1, xstart, nui, ui, rho, invrho, &
          eps6, eps8, eps10, combi, alp3, alp4, &
          n_vertices, n_facets, vertices_coor, facets_coor, n_extrpts, extrpts_coor 

!
end module commvars
