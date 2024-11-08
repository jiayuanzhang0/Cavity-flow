module mod_parameter
implicit none
 
! parameter
real*8, parameter   ::  Re = 1000d0  
 
! space
integer, parameter  ::  nx = 200
integer, parameter  ::  ny = 200
real*8,  parameter  ::  lx = 1d0
real*8,  parameter  ::  ly = 1d0

! time
real*8,  parameter  ::  dt = 1d-3
integer, parameter  ::  step_e = 100000     ! nsteps in total
integer, parameter  ::  inst_stride = 1000
integer, parameter  ::  resumable_stride = 20

! grid
real*8, dimension(0:nx)   ::  x
real*8, dimension(0:ny)   ::  y
real*8, dimension(1:nx)   ::  xc
real*8, dimension(1:ny)   ::  yc

! BC
real*8, dimension(0:nx  )   ::  u_n, u_s
real*8, dimension(0:ny  )   ::  v_w, v_e
real*8, dimension(0:nx+1)   ::  v_n, v_s
real*8, dimension(0:ny+1)   ::  u_w, u_e
 
! iteration
integer, parameter  ::  iter_max = 20
real*8,  parameter  ::  omg = 1.0d0   ! relaxation factor
real*8,  parameter  ::  tol = 1d-6

! dependent  
real*8,  parameter  ::  dx = lx/real(nx) 
real*8,  parameter  ::  dy = ly/real(ny)

! str
character(len=7)   ::  inst_dir = './inst/' 
character(len=12)  ::  resumable_dir = './resumable/' 
character(len=23)  ::  latestResumable_path = './resumable.Unformatted' 

! parallel
integer, parameter  :: num_threads = 4
 
contains
! grid
!------------------------------------------------------------------------------
subroutine get_grid(nx, ny, x, y, xc, yc)
implicit none

! in
integer, intent(in)  ::  nx
integer, intent(in)  ::  ny
! intermediate
integer ::  i, j
! out
real*8, dimension(0:nx), intent(out)    ::  x
real*8, dimension(0:ny), intent(out)    ::  y
real*8, dimension(1:nx), intent(out)    ::  xc
real*8, dimension(1:ny), intent(out)    ::  yc

do i = 0,nx
    x(i) = real(i)*dx
end do 
do j = 0,ny
    y(j) = real(j)*dy
end do 

do i = 1,nx
    xc(i) = (real(i)-0.5d0)*dx
end do 
do j = 1,ny
    yc(j) = (real(j)-0.5d0)*dy
end do
   
endsubroutine

! BC
!------------------------------------------------------------------------------
subroutine get_BC_mean(nx, ny, u_n, u_s, v_w, v_e, v_n, v_s, u_w, u_e)
implicit none

integer, intent(in)  ::  nx
integer, intent(in)  ::  ny

real*8, dimension(0:nx), intent(out)    ::  u_n, u_s
real*8, dimension(0:ny), intent(out)    ::  v_w, v_e
real*8, dimension(1:nx), intent(out)    ::  v_n, v_s
real*8, dimension(1:ny), intent(out)    ::  u_w, u_e
 
u_w = 0d0
u_e = 0d0
u_s = 0d0
u_n = 1d0

v_w = 0d0
v_e = 0d0
v_s = 0d0
v_n = 0d0
 
endsubroutine

endmodule 

