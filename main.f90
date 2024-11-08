! 2D Cavity Flow
! Created by Zhang Jia-Yuan on 2021.6.13 
! Acknowledge Gao Zhen-Yuan for debugging
! V1.4.10
! Projection method
! Time advance: RK2
! Pressure Poisson: SOR

program main
use mod_parameter
use omp_lib
implicit none

! declaration
!==============================================================================
! scalar
real*8  ::  t  
 
! array
real*8, dimension(0:nx  , 0:ny+1)   ::  u, u1, u_star, u_nxt, cfl_u  
real*8, dimension(0:nx+1, 0:ny  )   ::  v, v1, v_star, v_nxt, cfl_v
real*8, dimension(1:nx  , 1:ny  )   ::  uc, vc 
real*8, dimension(1:nx  , 1:ny  )   ::  p_Poisson_rhs
real*8, dimension(0:nx+1, 0:ny+1)   ::  p, p_new  

! iteration
real*8  ::  err

! intermediate
real*8  ::  p_temp  
integer ::  i, j, step, step_s, iter 

! post
real*8  ::  cfl, cfl_max, GAMA

! IO
logical alive

! parallel
integer :: proc_wts, proc_wte, proc_wt

! initialization
!==============================================================================
CALL omp_set_num_threads( num_threads )

! get velocity BC: u_n, u_s, v_w, v_e, v_n, v_s, u_w, u_e
! nx, ny cannot be reduced, because a subroutine cannot use the module it is in
CALL get_BC_mean(nx, ny, u_n, u_s, v_w, v_e, v_n, v_s, u_w, u_e)

! mkdir at beginning, so that the program will stop if dir is renamed
inquire(directory=inst_dir, exist=alive)
if (.not. alive) CALL system('mkdir '//inst_dir )
inquire(directory=resumable_dir, exist=alive)
if (.not. alive) CALL system('mkdir '//resumable_dir )

! initialization of fields
u      = 0d0    
v      = 0d0
u1     = 0d0
v1     = 0d0    
u_star = 0d0    ! *
v_star = 0d0    ! *
u_nxt  = 0d0    ! *
v_nxt  = 0d0    ! *
cfl_u  = 0d0
cfl_v  = 0d0
p      = 0d0    
p_new  = 0d0      

! initialization according to resumable file
inquire(file=latestResumable_path, exist=alive)
if(.not. alive) then
    write(*,*) 'Resumable file does not exist.'
    write(*,*) 'Create new case...'
    
    ! get grid info: x, y, xc, yc
    CALL get_grid(nx, ny, x, y, xc, yc)

    ! time
    t = 0d0
    step_s = 1

    ! set IC in bulk
    u = 0d0
    v = 0d0 
    p = 0d0 
    
    !! set BC
    ! set normal v
    CALL keepBC_vNor(u, v)
    ! complement tangential v at ghost region according to IC in the bulk
    CALL keepBC_vTan(u, v)    
    
    ! save initial field
    CALL cal_uc(u, v, uc, vc)
    CALL write_inst(step, t, uc, vc, p)
    
else
    open(11, file = latestResumable_path, access = 'stream', form = 'Unformatted') 
        write(*,*) 'Resumable file exists.'
        write(*,*) 'Read Resumable...'
        read(11) step_s, t, xc, yc, u, v, p
    close(11)
    step_s = step_s+1

end if

! initialize p_new for SOR
p_new = p
 
! calculation
!==============================================================================
write(*,*) 'Calculation starts!'
write(*,"('step_e = ', I9)")  step_e 
write(*,"('dt     = ', F40.20)")  dt 
 
CALL system_clock(proc_wts)     ! omp_get_wtime()
do step = step_s, step_e
    !write(*,"('step = ', I4)")  step

    ! substep 1: predetermine
    !--------------------------------------------------------------------------
    ! Eular
    !CALL predetermine(u, v, u_star, v_star)
    
    ! RK2
    CALL predetermine(u,v,u1,v1              )
    CALL keepBC_vTan(     u1,v1              )
    CALL keepBC_vNor(     u1,v1              )
    CALL predetermine(    u1,v1,u_star,v_star)
    ! CALL keepBC_vTan(           u_star,v_star)
    ! CALL keepBC_vNor(     u1,v1              )
    ! no need to keepBC for v, since the ghost region of v_star is not used
    ! keepBC_vNor for generalization: transient BC of vNor
 
    u_star = (u+u_star)/2d0
    v_star = (v+v_star)/2d0
 
    ! substep 2: pressure Poison 
    !--------------------------------------------------------------------------   
    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j)
    do j=1,ny         
    do i=1,nx
        p_Poisson_rhs(i,j) = ( ( u_star(i,j)-u_star(i-1,j) )/dx &       
                              +( v_star(i,j)-v_star(i,j-1) )/dy )/dt                        
    enddo
    enddo
    !$OMP END PARALLEL DO     
    
    do iter = 1,iter_max
        
        !! relaxation based on Jacobi
        !!keep Jacobi for parallel 
        !do j=1,ny
        !do i=1,nx
        !    p_temp = (p(i+1,j)+p(i-1,j))/(dx**2d0) + (p(i,j+1)+p(i,j-1))/(dy**2d0)
        !    p_new(i,j) = (1d0-omg)*p(i,j) & 
        !                + omg*(p_Poisson_rhs(i,j)-p_temp)/(-2d0/(dx**2d0)-2d0/(dy**2d0))
        !enddo
        !enddo
        
        ! relaxation based on Gauss-Seidel
        do j=1,ny
        do i=1,nx
            p_temp = (p_new(i+1,j)+p_new(i-1,j))/(dx**2d0) &
                    +(p_new(i,j+1)+p_new(i,j-1))/(dy**2d0)
            p_new(i,j) = (1d0-omg)*p_new(i,j) & 
                        +omg*(p_Poisson_rhs(i,j)-p_temp)/(-2d0/(dx**2d0)-2d0/(dy**2d0))
        enddo
        enddo 

        ! keep BC    
        CALL keepBC_p(p_new)
        
        ! cal err before update
        !err = sum(abs(p_new-p))/size(p) 
        err = maxval(abs(p_new-p))
        
        ! update p before exit
        p = p_new                       
        if(err < tol) exit  
    enddo
    
    ! iter correction
    if(iter==iter_max+1) then
        ! write(*,*) 'p not converge'
        iter=iter-1
    endif
    
    ! substep 3: final velocity 
    !--------------------------------------------------------------------------
    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j)
    do j=1,ny
    do i=0+1,nx-1 
        u_nxt(i,j) = u_star(i,j)-dt*(p(i+1,j)-p(i,j))/dx
    enddo
    enddo
    !$OMP END PARALLEL DO  
    
    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j)
    do j=0+1,ny-1
    do i=1,nx 
        v_nxt(i,j) = v_star(i,j)-dt*(p(i,j+1)-p(i,j))/dy
    enddo
    enddo   
    !$OMP END PARALLEL DO  
    
    ! update
    ! No need to update predetermine steps. Notice that RK is a single step method
    u = u_nxt
    v = v_nxt
    t = t+dt
    
    ! keep BC
    CALL keepBC_vTan(u, v)
    CALL keepBC_vNor(u, v)
    
    ! post
    !--------------------------------------------------------------------------
    ! cal uc first
    ! uc is used for integration in post processing
    CALL cal_uc(u, v, uc, vc)

    ! save inst
    if(mod(step, inst_stride).eq.0) then
        CALL write_inst(step, t, uc, vc, p)
    endif
    
    ! save resumable
    if(mod(step, resumable_stride).eq.0) then
        CALL write_resumable(step, t, u, v, p)
    endif
    
    !! save output
    ! cfl
    cfl_u = dabs(u)*dt/dx 
    cfl_v = dabs(v)*dt/dy 
    cfl_max = max( cfl_max, maxval(cfl_u), maxval(cfl_v) )
    ! GAMA
    CALL cal_GAMA(uc, vc, GAMA)
    ! save
    CALL write_output(step, iter, err, cfl_max, GAMA, u, v, p)
      
enddo

CALL system_clock(proc_wte) 
proc_wt  = proc_wte-proc_wts
write(*,*) 'proc_wt=', proc_wt
     
endprogram
 
 
!! subroutine
! calculation
!==============================================================================
! predetermine
!------------------------------------------------------------------------------
! the ghost region will not be calculated
subroutine predetermine(u, v, u_nxt, v_nxt)
use mod_parameter
use omp_lib
implicit none

! input
real*8, dimension(0:nx  , 0:ny+1), intent(in)   ::  u 
real*8, dimension(0:nx+1, 0:ny  ), intent(in)   ::  v 

! intermediate
real*8  ::  conv, visc
integer ::  i, j 
real*8, dimension(0:nx  , 0:ny  )   ::  un, vn

! output
real*8, dimension(0:nx  , 0:ny+1), intent(out)  ::  u_nxt
real*8, dimension(0:nx+1, 0:ny  ), intent(out)  ::  v_nxt

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j)
do j=0,ny
do i=0,nx
    un(i,j) = (u(i,j)+u(i,j+1))/2d0
    vn(i,j) = (v(i,j)+v(i+1,j))/2d0
enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j,conv,visc)
do j=1,ny
do i=0+1,nx-1 
    conv = ( u(i+1,j)**2d0-u(i-1,j)**2d0 )/(2d0*dx) &
          +( un(i,j)*vn(i,j) - un(i,j-1)*vn(i,j-1) )/dy
    visc = ( ( u(i+1,j)-2d0*u(i,j)+u(i-1,j) )/(dx**2d0) &       
            +( u(i,j+1)-2d0*u(i,j)+u(i,j-1) )/(dy**2d0) )/Re
    u_nxt(i,j) = u(i,j)+(visc-conv)*dt
enddo
enddo
!$OMP END PARALLEL DO 

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j,conv,visc)
do j=0+1,ny-1
do i=1,nx  
    conv = ( v(i,j+1)**2d0-v(i,j-1)**2d0 )/(2d0*dy) &
          +( un(i,j)*vn(i,j) - un(i-1,j)*vn(i-1,j) )/dy    
    visc = ( ( v(i+1,j)-2d0*v(i,j)+v(i-1,j) )/(dx**2d0) &       
            +( v(i,j+1)-2d0*v(i,j)+v(i,j-1) )/(dy**2d0) )/Re 
    v_nxt(i,j) = v(i,j)+(visc-conv)*dt
enddo
enddo
!$OMP END PARALLEL DO      
 
endsubroutine

! keep B.C. for tangential v
!------------------------------------------------------------------------------
subroutine keepBC_vTan(u, v)
use mod_parameter
implicit none

real*8, dimension(0:nx  , 0:ny+1),  intent(inout)   ::  u
real*8, dimension(0:nx+1, 0:ny  ),  intent(inout)   ::  v
 
u(:,   0) = -u(:, 1)+2d0*u_s 
u(:,ny+1) = -u(:,ny)+2d0*u_n   
v(0   ,:) = -v(1 ,:)+2d0*v_w 
v(nx+1,:) = -v(nx,:)+2d0*v_e 
 
endsubroutine 

! keep B.C. for normal v
!------------------------------------------------------------------------------
subroutine keepBC_vNor(u, v)
use mod_parameter
implicit none

real*8, dimension(0:nx  , 0:ny+1),  intent(inout)   ::  u
real*8, dimension(0:nx+1, 0:ny  ),  intent(inout)   ::  v
 
u(0 ,: ) = u_w
u(nx,: ) = u_e
v(: ,0 ) = v_s
v(: ,ny) = v_n
 
endsubroutine 
 
! keepBC_P
!------------------------------------------------------------------------------
subroutine keepBC_p(p)
use mod_parameter
implicit none

real*8, dimension(0:nx+1, 0:ny+1),  intent(inout)   ::  p
 
p(0   ,:) = p(1 ,:)
p(nx+1,:) = p(nx,:)
p(:,   0) = p(:, 1)
p(:,ny+1) = p(:,ny)
 
endsubroutine  

! cal_uc
!------------------------------------------------------------------------------
subroutine cal_uc(u, v, uc, vc)
use mod_parameter
implicit none

! input
real*8, dimension(0:nx  , 0:ny+1),  intent(in)   ::  u
real*8, dimension(0:nx+1, 0:ny  ),  intent(in)   ::  v

! intermediate
integer :: i,j

! output
real*8, dimension(1:nx  , 1:ny  ),  intent(out)  ::  uc, vc 

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j)
do j=1,ny
do i=1,nx     
    uc(i,j) = (u(i,j)+u(i-1,j))/2d0
    vc(i,j) = (v(i,j)+v(i,j-1))/2d0
enddo
enddo
!$OMP END PARALLEL DO 

endsubroutine 

! GAMA
!------------------------------------------------------------------------------
subroutine cal_GAMA(uc, vc, GAMA)
use mod_parameter
implicit none

! input
real*8, dimension(1:nx  , 1:ny  ),  intent(in)   ::  uc, vc

! intermediate
integer :: i, j

! output
real*8,  intent(out)   ::  GAMA

GAMA = 0d0

do j=1,ny
do i=1,nx
    GAMA = GAMA + ( uc(i,j)*(0.5d0*ly-yc(j)) + vc(i,j)*(xc(i)-0.5d0*lx) ) *dx*dy
enddo
enddo

endsubroutine 

! IO
!==============================================================================
! write inst
!------------------------------------------------------------------------------
subroutine write_inst(step, t, uc, vc, p)
use mod_parameter
implicit none

! input scalar
integer,    intent(in)  ::  step
real*8,     intent(in)  ::  t
! input array
real*8, dimension(1:nx  , 1:ny  ),  intent(in)   ::  uc, vc 
real*8, dimension(0:nx+1, 0:ny+1),  intent(in)   ::  p

! intermediate
integer :: i, j
character(len=100)  :: inst_path
character(len=100)  :: step_str

write(step_str, '(I9.9)') step   
inst_path = inst_dir//'inst_'//trim(adjustl(step_str))//'.Unformatted'

open(11, file = inst_path , access = 'stream' , form = 'Unformatted') 
    write(11) step, nx, ny, t, dt, xc, yc, uc, vc, p(1:nx, 1:ny)   
close(11) 
    
! formatted form; used for data check
!open(unit=11, file=inst_path, form='formatted')   
!    write(unit=11, fmt='(I32)') step, nx, ny
!    write(unit=11, fmt='(F32.16)') t, dt       
!    write(unit=11, fmt='(F32.16)') xc, yc 
!    write(unit=11, fmt='(F32.16)') uc, vc
!    write(unit=11, fmt='(F32.16)') p(1:nx, 1:ny)
!close(unit=11)

endsubroutine
 
! write resumable
!------------------------------------------------------------------------------
subroutine write_resumable(step, t, u, v, p)
use mod_parameter
implicit none

! input scalar
integer,    intent(in)  ::  step
real*8,     intent(in)  ::  t
! input array
real*8, dimension(0:nx  , 0:ny+1),  intent(in)   ::  u
real*8, dimension(0:nx+1, 0:ny  ),  intent(in)   ::  v
real*8, dimension(0:nx+1, 0:ny+1),  intent(in)   ::  p

! intermediate
character(len=100)  :: resumable_path
character(len=100)  :: step_str 

write(step_str, '(I9.9)') step   
resumable_path = resumable_dir//'resumable_'//trim(adjustl(step_str))//'.Unformatted'
 
open(11, file=resumable_path, access='stream', form='Unformatted') 
    write(11) step, t, xc, yc, u, v, p
close(11) 

CALL system('cp '//resumable_path//' '//latestResumable_path )
    
endsubroutine

! write output
!------------------------------------------------------------------------------
subroutine write_output(step, iter, err, cfl_max, GAMA, u, v, p)
use mod_parameter
implicit none

! input scalar
integer,    intent(in)  ::  step, iter
real*8,     intent(in)  ::  cfl_max, err, GAMA
! input array
real*8, dimension(0:nx  , 0:ny+1),  intent(in)   ::  u
real*8, dimension(0:nx+1, 0:ny  ),  intent(in)   ::  v
real*8, dimension(0:nx+1, 0:ny+1),  intent(in)   ::  p

! intermediate
real*8  ::  p_mean, u_max, u_min, v_max, v_min

p_mean = sum(p)/size(p)
u_max = maxval(u)
u_min = minval(u)
v_max = maxval(v)
v_min = minval(v)

open(11, file = 'output.txt', action='write', status='unknown', access='append') 
    ! NEVER forget to modify the number of output figures!
    write(11, "('step,', I, 4X, ',iter,', I, 4X, 8(',',A,',', F22.16, 4X))") &
                 step,            iter, &
                'p_mean'    ,p_mean ,&
                'err'       ,err    ,&
                'cfl_max'   ,cfl_max,&
                'u_max'     ,u_max  ,&
                'u_min'     ,u_min  ,&
                'v_max'     ,v_max  ,&
                'v_min'     ,v_min  ,&
                'GAMA'      ,GAMA
close(11) 

endsubroutine
