!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module scalar_flux

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : decomp_2d_abort
  use decomp_2d, only : xsize, ysize, zsize, ystart, zstart, &
       transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
  use variables, only : nx, ny, nz
  use param, only : dx, dy, dz, dx2, dy2, dz2, ifirstder, isecondder
  use derivX, only : alsaix, asix, bsix, csix, dsix
  use derivY, only : alsajy, asjy, bsjy, csjy, dsjy
  use derivZ, only : alsakz, askz, bskz, cskz, dskz
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

  implicit none
  private

  ! One unique face after each node, in its full directional pencil. These
  ! are flux densities: multiply by the actual face area for scalar mass.
  ! X face 1 is the physical inlet at x=0, and x face nx is the outlet.
  ! Prescribed node 1 is external data. Active node 2 owns width 3*dx/2,
  ! interiors own dx, and the outlet owns dx/2, covering the full domain.
  type, public :: scalar_faces
     real(mytype), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)
  end type scalar_faces

  type :: line_factors
     integer :: n = 0
     logical :: cyclic = .false.
     real(mytype) :: corner_ratio = 0._mytype, correction_denominator = 1._mytype
     real(mytype), allocatable :: upper(:), multiplier(:), inverse_diagonal(:), correction(:)
  end type line_factors

  type(line_factors), save :: convection(3), diffusion(3)
  real(mytype), save :: diffusion_weights(4,3)
  logical, save :: initialized = .false.
  real(mytype), allocatable, save :: work_y(:,:,:), velocity_y(:,:,:)
  real(mytype), allocatable, save :: work_z(:,:,:), velocity_z(:,:,:), work_x(:,:,:)

  public :: faces_allocate, faces_zero, faces_copy, faces_deallocate
  public :: scalar_raw_velocity, scalar_high_flux, scalar_divergence
  public :: scalar_cell_volume, scalar_cell_width

  ! Reuse the solver's existing generic tridiagonal factorization.
  interface
     subroutine prepare(b,c,f,s,w,n)
       import :: mytype
       integer, intent(in) :: n
       real(mytype), intent(in) :: b(n), c(n), f(n)
       real(mytype), intent(out) :: s(n), w(n)
     end subroutine prepare
  end interface

contains

  subroutine faces_allocate(f)
    type(scalar_faces), intent(inout) :: f
    if (.not.allocated(f%x)) allocate(f%x(xsize(1),xsize(2),xsize(3)))
    if (.not.allocated(f%y)) allocate(f%y(ysize(1),ysize(2),ysize(3)))
    if (.not.allocated(f%z)) allocate(f%z(zsize(1),zsize(2),zsize(3)))
  end subroutine faces_allocate

  subroutine faces_deallocate(f)
    type(scalar_faces), intent(inout) :: f
    if (allocated(f%x)) deallocate(f%x)
    if (allocated(f%y)) deallocate(f%y)
    if (allocated(f%z)) deallocate(f%z)
  end subroutine faces_deallocate

  subroutine faces_zero(f)
    type(scalar_faces), intent(inout) :: f
    call faces_allocate(f)
    f%x = 0._mytype
    f%y = 0._mytype
    f%z = 0._mytype
  end subroutine faces_zero

  subroutine faces_copy(a,b)
    type(scalar_faces), intent(in) :: a
    type(scalar_faces), intent(inout) :: b
    call faces_allocate(b)
    b%x = a%x
    b%y = a%y
    b%z = a%z
  end subroutine faces_copy

  function scalar_cell_width(global_i) result(width)
    integer, intent(in) :: global_i
    real(mytype) :: width
    width = dx
    if (global_i == 1) width = 0._mytype
    if (global_i == 2) width = 1.5_mytype*dx
    if (global_i == nx) width = 0.5_mytype*dx
  end function scalar_cell_width

  function scalar_cell_volume(global_i) result(volume)
    integer, intent(in) :: global_i
    real(mytype) :: volume
    volume = scalar_cell_width(global_i)*dy*dz
  end function scalar_cell_volume

  subroutine initialize_flux()
    real(mytype) :: alpha(3)
    integer :: d, lengths(3)

    if (initialized) return
    if (ifirstder /= 4 .or. isecondder /= 5) &
         call decomp_2d_abort(1,'Conservative scalar flux requires ifirstder=4, isecondder=5')
    if (nx < 10 .or. ny < 8 .or. nz < 8) &
         call decomp_2d_abort(1,'Conservative scalar flux requires nx>=10, ny/nz>=8')
    if (min(dx,dy,dz) <= 0._mytype) &
         call decomp_2d_abort(1,'Invalid conservative scalar mesh spacing')

    ! schemes() has already generated these coefficients from the current
    ! configuration. Sharing them avoids a second iSVV coefficient generator.
    alpha = [alsaix, alsajy, alsakz]
    diffusion_weights(:,1) = [asix,bsix,csix,dsix]*dx2
    diffusion_weights(:,2) = [asjy,bsjy,csjy,dsjy]*dy2
    diffusion_weights(:,3) = [askz,bskz,cskz,dskz]*dz2
    if (.not.all(ieee_is_finite(alpha)) .or. .not.all(ieee_is_finite(diffusion_weights))) &
         call decomp_2d_abort(1,'Nonfinite conservative scalar iSVV coefficients')
    ! This first implementation accepts the nonsingular diagonally dominant
    ! compact family only. Failure cannot silently choose another operator.
    if (any(abs(alpha) >= 0.5_mytype)) &
         call decomp_2d_abort(1,'Conservative scalar iSVV compact matrix is not diagonally dominant')
    lengths = [nx,ny,nz]
    do d=1,3
       call factor_line(convection(d),lengths(d),1._mytype/3._mytype,d/=1,.false.)
       call factor_line(diffusion(d),lengths(d),alpha(d),d/=1,.true.)
    enddo

    allocate(work_y(ysize(1),ysize(2),ysize(3)),velocity_y(ysize(1),ysize(2),ysize(3)))
    allocate(work_z(zsize(1),zsize(2),zsize(3)),velocity_z(zsize(1),zsize(2),zsize(3)))
    allocate(work_x(xsize(1),xsize(2),xsize(3)))
    initialized = .true.
  end subroutine initialize_flux

  subroutine factor_line(f,n,alpha,periodic,is_diffusion)
    type(line_factors), intent(inout) :: f
    integer, intent(in) :: n
    real(mytype), intent(in) :: alpha
    logical, intent(in) :: periodic, is_diffusion
    real(mytype) :: lower(n), diagonal(n), upper(n), rhs(n), gamma
    integer :: i,m

    m=n
    if (.not.periodic) m=n-1
    f%n=m
    f%cyclic=periodic
    allocate(f%upper(m),f%multiplier(m),f%inverse_diagonal(m),f%correction(m))
    lower=0._mytype
    upper=0._mytype
    diagonal=1._mytype
    do i=1,m
       if (.not.periodic) then
          if (is_diffusion) then
             if (i<4 .or. i>n-4) cycle
          else
             if (i==1 .or. i==m) cycle
          endif
       endif
       if (i>1) lower(i-1)=alpha
       if (i<m) upper(i)=alpha
    enddo
    if (periodic) then
       gamma=-1._mytype
       diagonal(1)=diagonal(1)-gamma
       diagonal(m)=diagonal(m)-alpha*alpha/gamma
       f%corner_ratio=alpha/gamma
    endif
    f%upper=upper(1:m)
    f%multiplier=0._mytype
    call prepare(lower(1:m),diagonal(1:m),upper(1:m),f%multiplier,f%inverse_diagonal,m)
    f%correction=0._mytype
    if (periodic) then
       rhs=0._mytype
       rhs(1)=gamma
       rhs(m)=alpha
       call solve_tridiagonal(f,rhs(1:m),lower(1:m))
       f%correction=lower(1:m)
       f%correction_denominator=1._mytype+f%correction(1)+f%corner_ratio*f%correction(m)
       if (abs(f%correction_denominator)<epsilon(1._mytype)) &
            call decomp_2d_abort(1,'Singular conservative scalar cyclic compact matrix')
    endif
  end subroutine factor_line

  subroutine solve_tridiagonal(f,rhs,solution)
    type(line_factors), intent(in) :: f
    real(mytype), intent(in) :: rhs(:)
    real(mytype), intent(out) :: solution(:)
    integer :: i
    solution(1)=rhs(1)
    do i=2,f%n
       solution(i)=rhs(i)-f%multiplier(i)*solution(i-1)
    enddo
    solution(f%n)=solution(f%n)*f%inverse_diagonal(f%n)
    do i=f%n-1,1,-1
       solution(i)=(solution(i)-f%upper(i)*solution(i+1))*f%inverse_diagonal(i)
    enddo
  end subroutine solve_tridiagonal

  subroutine solve_line(f,rhs,solution)
    type(line_factors), intent(in) :: f
    real(mytype), intent(in) :: rhs(:)
    real(mytype), intent(out) :: solution(:)
    real(mytype) :: correction
    call solve_tridiagonal(f,rhs,solution)
    if (f%cyclic) then
       correction=(solution(1)+f%corner_ratio*solution(f%n))/f%correction_denominator
       solution(1:f%n)=solution(1:f%n)-correction*f%correction
    endif
  end subroutine solve_line

  pure integer function wrap(i,n)
    integer, intent(in) :: i,n
    wrap=modulo(i-1,n)+1
  end function wrap

  subroutine reconstruct_line(q,d,face)
    real(mytype), intent(in) :: q(:)
    integer, intent(in) :: d
    real(mytype), intent(out) :: face(:)
    real(mytype) :: rhs(size(q))
    integer :: i,n,m
    n=size(q)
    m=convection(d)%n
    if (d==1) then
       rhs(1)=0.5_mytype*(q(1)+q(2))
       rhs(m)=0.5_mytype*(q(n-1)+q(n))
       do i=2,m-1
          rhs(i)=(q(i-1)+29._mytype*q(i)+29._mytype*q(i+1)+q(i+2))/36._mytype
       enddo
    else
       do i=1,n
          rhs(i)=(q(wrap(i-1,n))+29._mytype*q(i)+29._mytype*q(wrap(i+1,n)) &
               +q(wrap(i+2,n)))/36._mytype
       enddo
    endif
    call solve_line(convection(d),rhs(1:m),face(1:m))
    if (d==1) face(n)=q(n)
  end subroutine reconstruct_line

  subroutine diffusion_line(q,d,h,face)
    real(mytype), intent(in) :: q(:),h
    integer, intent(in) :: d
    real(mytype), intent(out) :: face(:)
    real(mytype) :: rhs(size(q))
    integer :: i,n,m,r,s,il,ir
    n=size(q)
    m=diffusion(d)%n
    rhs=0._mytype
    do i=1,m
       if (d==1) then
          if (i<4 .or. i>n-4) then
             rhs(i)=(q(i+1)-q(i))/h
             cycle
          endif
       endif
       do r=1,4
          do s=0,r-1
             ir=i+1+s
             il=i-s
             if (d/=1) then
                ir=wrap(ir,n)
                il=wrap(il,n)
             endif
             rhs(i)=rhs(i)+diffusion_weights(r,d)*(q(ir)-q(il))/h
          enddo
       enddo
    enddo
    call solve_line(diffusion(d),rhs(1:m),face(1:m))
    if (d==1) face(n)=(3._mytype*q(n)-4._mytype*q(n-1)+q(n-2))/(2._mytype*h)
  end subroutine diffusion_line

  subroutine high_line(q,u,mapped,kappa,d,h,face)
    real(mytype), intent(in) :: q(:),u(:),mapped(:),kappa,h
    integer, intent(in) :: d
    real(mytype), intent(out) :: face(:)
    real(mytype) :: rq(size(q)),ru(size(q)),ruq(size(q)),gradient(size(q))
    call reconstruct_line(q,d,rq)
    call reconstruct_line(u,d,ru)
    call reconstruct_line(u*q,d,ruq)
    call diffusion_line(q,d,h,gradient)
    ! The explicit auxiliary-velocity defect vanishes without mapping. It
    ! preserves constants but is a distinct candidate when mapping is nonzero.
    face=ruq+(mapped-ru)*rq-kappa*gradient
  end subroutine high_line

  subroutine scalar_raw_velocity(ux,uy,uz,flow)
    real(mytype), intent(in) :: ux(xsize(1),xsize(2),xsize(3))
    real(mytype), intent(in) :: uy(xsize(1),xsize(2),xsize(3))
    real(mytype), intent(in) :: uz(xsize(1),xsize(2),xsize(3))
    type(scalar_faces), intent(inout) :: flow
    integer :: i,j,k
    call initialize_flux()
    call faces_allocate(flow)
    do k=1,xsize(3)
       do j=1,xsize(2)
          call reconstruct_line(ux(:,j,k),1,flow%x(:,j,k))
       enddo
    enddo
    flow%x(1,:,:)=ux(1,:,:)
    flow%x(nx,:,:)=ux(nx,:,:)
    call transpose_x_to_y(uy,velocity_y)
    do k=1,ysize(3)
       do i=1,ysize(1)
          if (i+ystart(1)-1==1) then
             flow%y(i,:,k)=0._mytype
          else
             call reconstruct_line(velocity_y(i,:,k),2,flow%y(i,:,k))
          endif
       enddo
    enddo
    call transpose_x_to_y(uz,velocity_y)
    call transpose_y_to_z(velocity_y,velocity_z)
    do j=1,zsize(2)
       do i=1,zsize(1)
          if (i+zstart(1)-1==1) then
             flow%z(i,j,:)=0._mytype
          else
             call reconstruct_line(velocity_z(i,j,:),3,flow%z(i,j,:))
          endif
       enddo
    enddo
  end subroutine scalar_raw_velocity

  subroutine scalar_high_flux(q,ux,uy,uz,mapped,kappa,flow)
    real(mytype), intent(in) :: q(xsize(1),xsize(2),xsize(3))
    real(mytype), intent(in) :: ux(xsize(1),xsize(2),xsize(3))
    real(mytype), intent(in) :: uy(xsize(1),xsize(2),xsize(3))
    real(mytype), intent(in) :: uz(xsize(1),xsize(2),xsize(3)),kappa
    type(scalar_faces), intent(in) :: mapped
    type(scalar_faces), intent(inout) :: flow
    integer :: i,j,k
    call initialize_flux()
    call faces_allocate(flow)
    do k=1,xsize(3)
       do j=1,xsize(2)
          call high_line(q(:,j,k),ux(:,j,k),mapped%x(:,j,k),kappa,1,dx,flow%x(:,j,k))
       enddo
    enddo
    call transpose_x_to_y(q,work_y)
    call transpose_x_to_y(uy,velocity_y)
    do k=1,ysize(3)
       do i=1,ysize(1)
          if (i+ystart(1)-1==1) then
             flow%y(i,:,k)=0._mytype
          else
             call high_line(work_y(i,:,k),velocity_y(i,:,k),mapped%y(i,:,k),kappa,2,dy,flow%y(i,:,k))
          endif
       enddo
    enddo
    call transpose_y_to_z(work_y,work_z)
    call transpose_x_to_y(uz,velocity_y)
    call transpose_y_to_z(velocity_y,velocity_z)
    do j=1,zsize(2)
       do i=1,zsize(1)
          if (i+zstart(1)-1==1) then
             flow%z(i,j,:)=0._mytype
          else
             call high_line(work_z(i,j,:),velocity_z(i,j,:),mapped%z(i,j,:),kappa,3,dz,flow%z(i,j,:))
          endif
       enddo
    enddo
  end subroutine scalar_high_flux

  subroutine scalar_divergence(flow,div)
    type(scalar_faces), intent(in) :: flow
    real(mytype), intent(out) :: div(xsize(1),xsize(2),xsize(3))
    integer :: i,j,k,jm,km
    call initialize_flux()
    do k=1,ysize(3)
       do j=1,ysize(2)
          jm=wrap(j-1,ny)
          work_y(:,j,k)=(flow%y(:,j,k)-flow%y(:,jm,k))/dy
       enddo
    enddo
    do k=1,zsize(3)
       km=wrap(k-1,nz)
       work_z(:,:,k)=(flow%z(:,:,k)-flow%z(:,:,km))/dz
    enddo
    call transpose_z_to_y(work_z,velocity_y)
    work_y=work_y+velocity_y
    call transpose_y_to_x(work_y,work_x)
    div=0._mytype
    do i=2,nx
       div(i,:,:)=(flow%x(i,:,:)-flow%x(i-1,:,:))/scalar_cell_width(i)+work_x(i,:,:)
    enddo
  end subroutine scalar_divergence

end module scalar_flux
