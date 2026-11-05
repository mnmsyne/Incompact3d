!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module scalar_fct

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : real_type, decomp_2d_abort, nrank
  use decomp_2d, only : xsize, ysize, zsize, ystart, zstart, &
       transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
  use mpi
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use variables, only : nx, ny, nz, numscalar, sc
  use param, only : dx, dy, dz, dt, xnu, itr, itime, ntime, adt, bdt, gdt, &
       scalar_scheme, scalar_fct_lower, scalar_fct_upper, sync_scal_needed
  use scalar_flux
  use scalar_flux_mapping, only : scalar_map_velocity

  implicit none
  private
  public :: scalar_fct_rhs, scalar_fct_advance

  ! Scalars and budgets remain in local pencils. Node 1 is prescribed inlet
  ! data, never an evolution unknown. The adjacent active volume includes
  ! its half-cell, so the mass budget covers the whole physical domain.
  logical, save :: initialized = .false.
  integer, save :: pending_stage = 0
  type(scalar_faces), save :: mapped, mass_integral, stage_flux, accepted_integral
  type(scalar_faces), allocatable, save :: high_integral(:), low_integral(:)
  real(mytype), allocatable, save :: old(:,:,:,:), low(:,:,:,:)
  real(mytype), allocatable, save :: div_x(:,:,:), plus_x(:,:,:), minus_x(:,:,:)
  real(mytype), allocatable, save :: ratio_plus(:,:,:), ratio_minus(:,:,:)
  real(mytype), allocatable, save :: q_y(:,:,:), q_z(:,:,:), work_y(:,:,:), work_z(:,:,:), work_x(:,:,:)
  real(mytype), allocatable, save :: rp_y(:,:,:), rm_y(:,:,:), rp_z(:,:,:), rm_z(:,:,:)
  real(mytype), allocatable, save :: cfl_x(:,:,:)

contains

  logical function new_scalar(is)
    integer, intent(in) :: is
    new_scalar = scalar_scheme(is)=='div'.or.scalar_scheme(is)=='fct'
  end function new_scalar

  subroutine initialize_fct()
    integer :: is
    if (initialized) return
    call faces_allocate(mapped)
    call faces_allocate(mass_integral)
    call faces_allocate(stage_flux)
    call faces_allocate(accepted_integral)
    allocate(high_integral(numscalar),low_integral(numscalar))
    do is=1,numscalar
       if (.not.new_scalar(is)) cycle
       call faces_allocate(high_integral(is))
       call faces_allocate(low_integral(is))
    enddo
    allocate(old(xsize(1),xsize(2),xsize(3),numscalar),low(xsize(1),xsize(2),xsize(3),numscalar))
    allocate(div_x(xsize(1),xsize(2),xsize(3)),plus_x(xsize(1),xsize(2),xsize(3)))
    allocate(minus_x(xsize(1),xsize(2),xsize(3)),ratio_plus(xsize(1),xsize(2),xsize(3)))
    allocate(ratio_minus(xsize(1),xsize(2),xsize(3)),cfl_x(xsize(1),xsize(2),xsize(3)))
    allocate(work_x(xsize(1),xsize(2),xsize(3)))
    allocate(q_y(ysize(1),ysize(2),ysize(3)),q_z(zsize(1),zsize(2),zsize(3)))
    allocate(work_y(ysize(1),ysize(2),ysize(3)),work_z(zsize(1),zsize(2),zsize(3)))
    allocate(rp_y(ysize(1),ysize(2),ysize(3)),rm_y(ysize(1),ysize(2),ysize(3)))
    allocate(rp_z(zsize(1),zsize(2),zsize(3)),rm_z(zsize(1),zsize(2),zsize(3)))
    initialized = .true.
  end subroutine initialize_fct

  subroutine check_finite(q,label)
    real(mytype), intent(in) :: q(:,:,:)
    character(len=*), intent(in) :: label
    integer :: bad, global_bad, ierr
    bad=0
    if (.not.all(ieee_is_finite(q))) bad=1
    call MPI_ALLREDUCE(bad,global_bad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (ierr/=MPI_SUCCESS) call decomp_2d_abort(1,'Scalar flux MPI check failed')
    if (global_bad/=0) call decomp_2d_abort(1,label)
  end subroutine check_finite

  subroutine global_extrema(q,bounds)
    real(mytype), intent(in) :: q(:,:,:)
    real(mytype), intent(out) :: bounds(2)
    real(mytype) :: local(2), both(2)
    integer :: ierr
    local=[-minval(q(2:nx,:,:)),maxval(q(2:nx,:,:))]
    call MPI_ALLREDUCE(local,both,2,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (ierr/=MPI_SUCCESS) call decomp_2d_abort(1,'Scalar flux extrema reduction failed')
    bounds=[-both(1),both(2)]
  end subroutine global_extrema

  subroutine add_faces(target,source,weight)
    type(scalar_faces), intent(inout) :: target
    type(scalar_faces), intent(in) :: source
    real(mytype), intent(in) :: weight
    target%x=target%x+weight*source%x
    target%y=target%y+weight*source%y
    target%z=target%z+weight*source%z
  end subroutine add_faces

  ! Called from transeq:scalar while the actual stage input velocities still
  ! exist. No velocity is reconstructed from the later momentum/projection.
  subroutine scalar_fct_rhs(phi,ux,uy,uz,dphi)
    real(mytype), intent(in) :: phi(xsize(1),xsize(2),xsize(3),numscalar)
    real(mytype), intent(in) :: ux(xsize(1),xsize(2),xsize(3)),uy(xsize(1),xsize(2),xsize(3))
    real(mytype), intent(in) :: uz(xsize(1),xsize(2),xsize(3))
    real(mytype), intent(inout) :: dphi(xsize(1),xsize(2),xsize(3),ntime,numscalar)
    integer :: is
    real(mytype) :: weight,kappa,bounds(2),tol
    call initialize_fct()
    if (itr/=pending_stage+1) call decomp_2d_abort(1,'Invalid conservative scalar RK stage ordering')
    tol=128._mytype*epsilon(1._mytype)*max(1._mytype,abs(scalar_fct_lower),abs(scalar_fct_upper))
    if (itr==1) then
       old=phi
       call faces_zero(mass_integral)
       do is=1,numscalar
          if (.not.new_scalar(is)) cycle
          call check_finite(old(:,:,:,is),'Nonfinite accepted conservative scalar state')
          if (scalar_scheme(is)=='fct') then
             call global_extrema(old(:,:,:,is),bounds)
             if (bounds(1)<scalar_fct_lower-tol.or.bounds(2)>scalar_fct_upper+tol) &
                  call decomp_2d_abort(1,'Accepted scalar is outside the FCT interval')
          endif
          call faces_zero(high_integral(is))
          dphi(:,:,:,2,is)=0._mytype
       enddo
    endif
    call scalar_map_velocity(ux,uy,uz,mapped)
    weight=0._mytype
    if (itr==1) weight=dt/4._mytype
    if (itr==3) weight=3._mytype*dt/4._mytype
    call add_faces(mass_integral,mapped,weight)
    do is=1,numscalar
       if (.not.new_scalar(is)) cycle
       call check_finite(phi(:,:,:,is),'Nonfinite conservative scalar RK stage')
       kappa=xnu/sc(is)
       ! Only the current and past stages are needed. The complete baseline
       ! and its joint CFL/range are checked BEFORE the affected stage-3 RHS.
       if (itr==3) call make_low_candidate(is,kappa)
       call scalar_high_flux(phi(:,:,:,is),ux,uy,uz,mapped,kappa,stage_flux)
       ! First version deliberately uses a low-order physical outlet baseline,
       ! without extrapolating the scalar state. The inlet stays Dirichlet.
       ! Both candidates use these identical physical time-integrated fluxes.
       stage_flux%x(1,:,:)=mapped%x(1,:,:)*old(1,:,:,is) &
            -kappa*(old(2,:,:,is)-old(1,:,:,is))/dx
       stage_flux%x(nx,:,:)=mapped%x(nx,:,:)*old(nx,:,:,is)
       call scalar_divergence(stage_flux,div_x)
       dphi(:,:,:,1,is)=-div_x
       call check_finite(dphi(:,:,:,1,is),'Nonfinite conservative scalar RHS')
       call add_faces(high_integral(is),stage_flux,weight)
    enddo
    pending_stage=itr
  end subroutine scalar_fct_rhs

  subroutine make_low_candidate(is,kappa)
    integer, intent(in) :: is
    real(mytype), intent(in) :: kappa
    integer :: i,j,k,jp,jm,kp,km,gi,ierr
    real(mytype) :: m,delta,bounds(2),cfl,local_cfl,tol
    type(scalar_faces) :: f
    call faces_allocate(f)
    f%x=0._mytype; f%y=0._mytype; f%z=0._mytype
    cfl_x=0._mytype
    do k=1,xsize(3)
       do j=1,xsize(2)
          f%x(1,j,k)=mass_integral%x(1,j,k)*old(1,j,k,is) &
               -dt*kappa*(old(2,j,k,is)-old(1,j,k,is))/dx
          do i=2,nx-1
             m=mass_integral%x(i,j,k)
             f%x(i,j,k)=max(m,0._mytype)*old(i,j,k,is)+min(m,0._mytype)*old(i+1,j,k,is) &
                  -dt*kappa*(old(i+1,j,k,is)-old(i,j,k,is))/dx
          enddo
          f%x(nx,j,k)=mass_integral%x(nx,j,k)*old(nx,j,k,is)
          do i=2,nx
             delta=dt*kappa/dx
             if (i<nx) delta=2._mytype*delta
             cfl_x(i,j,k)=(max(mass_integral%x(i,j,k),0._mytype) &
                  +max(-mass_integral%x(i-1,j,k),0._mytype)+delta)/scalar_cell_width(i)
          enddo
       enddo
    enddo
    call transpose_x_to_y(old(:,:,:,is),q_y)
    call transpose_y_to_z(q_y,q_z)
    work_y=0._mytype; work_z=0._mytype
    do k=1,ysize(3)
       do j=1,ny
          jp=modulo(j,ny)+1; jm=modulo(j-2,ny)+1
          do i=1,ysize(1)
             gi=i+ystart(1)-1
             if (gi==1) cycle
             m=mass_integral%y(i,j,k)
             f%y(i,j,k)=max(m,0._mytype)*q_y(i,j,k)+min(m,0._mytype)*q_y(i,jp,k) &
                  -dt*kappa*(q_y(i,jp,k)-q_y(i,j,k))/dy
             work_y(i,j,k)=(max(m,0._mytype)+max(-mass_integral%y(i,jm,k),0._mytype) &
                  +2._mytype*dt*kappa/dy)/dy
          enddo
       enddo
    enddo
    do k=1,nz
       kp=modulo(k,nz)+1; km=modulo(k-2,nz)+1
       do j=1,zsize(2)
          do i=1,zsize(1)
             gi=i+zstart(1)-1
             if (gi==1) cycle
             m=mass_integral%z(i,j,k)
             f%z(i,j,k)=max(m,0._mytype)*q_z(i,j,k)+min(m,0._mytype)*q_z(i,j,kp) &
                  -dt*kappa*(q_z(i,j,kp)-q_z(i,j,k))/dz
             work_z(i,j,k)=(max(m,0._mytype)+max(-mass_integral%z(i,j,km),0._mytype) &
                  +2._mytype*dt*kappa/dz)/dz
          enddo
       enddo
    enddo
    call transpose_y_to_x(work_y,work_x)
    cfl_x=cfl_x+work_x
    call transpose_z_to_y(work_z,work_y)
    call transpose_y_to_x(work_y,work_x)
    cfl_x=cfl_x+work_x
    local_cfl=maxval(cfl_x)
    call MPI_ALLREDUCE(local_cfl,cfl,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (ierr/=MPI_SUCCESS) call decomp_2d_abort(1,'Scalar CFL reduction failed')
    call faces_copy(f,low_integral(is))
    call faces_deallocate(f)
    call scalar_divergence(low_integral(is),div_x)
    low(:,:,:,is)=old(:,:,:,is)-div_x
    call check_finite(low(:,:,:,is),'Nonfinite low-order scalar candidate')
    call global_extrema(low(:,:,:,is),bounds)
    if (nrank==0) write(*,*) 'SCALAR_BASE step,species,CFL,min,max=',itime,is,cfl,bounds
    if (scalar_scheme(is)=='fct') then
       tol=256._mytype*epsilon(1._mytype)*max(1._mytype,abs(scalar_fct_lower),abs(scalar_fct_upper))
       if (cfl>1._mytype+tol) call decomp_2d_abort(1,'Joint scalar advection/diffusion CFL exceeds one')
       if (bounds(1)<scalar_fct_lower-tol.or.bounds(2)>scalar_fct_upper+tol) &
            call decomp_2d_abort(1,'Baseline low-order complete scalar step is inadmissible')
    endif
    ! No predictor increment in this first production closure: thetaB=1,
    ! deltaI=0. The baseline itself has already been checked when FCT is used.
    if (nrank==0) write(*,*) 'SCALAR_BOUNDARY species,thetaB,correction_net,L1,max=', &
         is,1._mytype,0._mytype,0._mytype,0._mytype
  end subroutine make_low_candidate

  ! Accumulate one face's two opposite mass exchanges in a full pencil.
  subroutine add_margin(a,pl,ml,pr,mr)
    real(mytype), intent(in) :: a
    real(mytype), intent(inout) :: pl,ml,pr,mr
    pl=pl+max(-a,0._mytype); ml=ml+max(a,0._mytype)
    pr=pr+max(a,0._mytype); mr=mr+max(-a,0._mytype)
  end subroutine add_margin

  subroutine compute_ratios(anti,is)
    type(scalar_faces), intent(in) :: anti
    integer, intent(in) :: is
    integer :: i,j,k,jp,kp,gi
    real(mytype) :: area,volume
    plus_x=0._mytype; minus_x=0._mytype
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=2,nx-1
             call add_margin(anti%x(i,j,k)*dy*dz,plus_x(i,j,k),minus_x(i,j,k), &
                  plus_x(i+1,j,k),minus_x(i+1,j,k))
          enddo
       enddo
    enddo
    rp_y=0._mytype; rm_y=0._mytype
    do k=1,ysize(3)
       do j=1,ny
          jp=modulo(j,ny)+1
          do i=1,ysize(1)
             gi=i+ystart(1)-1
             if (gi==1) cycle
             area=scalar_cell_width(gi)*dz
             call add_margin(anti%y(i,j,k)*area,rp_y(i,j,k),rm_y(i,j,k),rp_y(i,jp,k),rm_y(i,jp,k))
          enddo
       enddo
    enddo
    call transpose_y_to_x(rp_y,work_x); plus_x=plus_x+work_x
    call transpose_y_to_x(rm_y,work_x); minus_x=minus_x+work_x
    rp_z=0._mytype; rm_z=0._mytype
    do k=1,nz
       kp=modulo(k,nz)+1
       do j=1,zsize(2)
          do i=1,zsize(1)
             gi=i+zstart(1)-1
             if (gi==1) cycle
             area=scalar_cell_width(gi)*dy
             call add_margin(anti%z(i,j,k)*area,rp_z(i,j,k),rm_z(i,j,k),rp_z(i,j,kp),rm_z(i,j,kp))
          enddo
       enddo
    enddo
    call transpose_z_to_y(rp_z,work_y)
    call transpose_y_to_x(work_y,work_x); plus_x=plus_x+work_x
    call transpose_z_to_y(rm_z,work_y)
    call transpose_y_to_x(work_y,work_x); minus_x=minus_x+work_x
    ratio_plus=1._mytype; ratio_minus=1._mytype
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=2,nx
             volume=scalar_cell_volume(i)
             if (plus_x(i,j,k)>0._mytype) ratio_plus(i,j,k)=min(1._mytype, &
                  max(0._mytype,(scalar_fct_upper-low(i,j,k,is))*volume)/plus_x(i,j,k))
             if (minus_x(i,j,k)>0._mytype) ratio_minus(i,j,k)=min(1._mytype, &
                  max(0._mytype,(low(i,j,k,is)-scalar_fct_lower)*volume)/minus_x(i,j,k))
          enddo
       enddo
    enddo
    call transpose_x_to_y(ratio_plus,rp_y); call transpose_x_to_y(ratio_minus,rm_y)
    call transpose_y_to_z(rp_y,rp_z); call transpose_y_to_z(rm_y,rm_z)
  end subroutine compute_ratios

  real(mytype) function face_theta(a,pl,ml,pr,mr) result(theta)
    real(mytype), intent(in) :: a,pl,ml,pr,mr
    if (a>=0._mytype) then
       theta=min(ml,pr)
    else
       theta=min(pl,mr)
    endif
  end function face_theta

  subroutine limit_internal(is,theta_min,removed_l1,removed_max,nlimited)
    integer, intent(in) :: is
    real(mytype), intent(out) :: theta_min,removed_l1,removed_max
    integer, intent(out) :: nlimited
    integer :: i,j,k,jp,kp,gi,ierr,count_local
    real(mytype) :: theta,area,local(3),global(3),a,boundary_error,global_error,scale,global_scale
    ! Reuse accepted_integral as the antidiffusive difference before limiting.
    accepted_integral%x=high_integral(is)%x-low_integral(is)%x
    accepted_integral%y=high_integral(is)%y-low_integral(is)%y
    accepted_integral%z=high_integral(is)%z-low_integral(is)%z
    boundary_error=max(maxval(abs(accepted_integral%x(1,:,:))),maxval(abs(accepted_integral%x(nx,:,:))))
    call MPI_ALLREDUCE(boundary_error,global_error,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (ierr/=MPI_SUCCESS) call decomp_2d_abort(1,'Scalar physical-flux reduction failed')
    scale=max(1._mytype,maxval(abs(low_integral(is)%x(1,:,:))),maxval(abs(low_integral(is)%x(nx,:,:))))
    call MPI_ALLREDUCE(scale,global_scale,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (ierr/=MPI_SUCCESS) call decomp_2d_abort(1,'Scalar physical-flux scale reduction failed')
    if (global_error>1024._mytype*epsilon(1._mytype)*global_scale) &
         call decomp_2d_abort(1,'High/low scalar physical boundary integrals disagree')
    accepted_integral%x(1,:,:)=0._mytype; accepted_integral%x(nx,:,:)=0._mytype
    call compute_ratios(accepted_integral,is)
    local=[1._mytype,0._mytype,0._mytype]; count_local=0
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=2,nx-1
             a=accepted_integral%x(i,j,k)
             theta=face_theta(a,ratio_plus(i,j,k),ratio_minus(i,j,k), &
                  ratio_plus(i+1,j,k),ratio_minus(i+1,j,k))
             call record(theta,a*dy*dz)
             accepted_integral%x(i,j,k)=theta*a
          enddo
       enddo
    enddo
    do k=1,ysize(3)
       do j=1,ny
          jp=modulo(j,ny)+1
          do i=1,ysize(1)
             gi=i+ystart(1)-1
             if (gi==1) cycle
             a=accepted_integral%y(i,j,k)
             theta=face_theta(a,rp_y(i,j,k),rm_y(i,j,k),rp_y(i,jp,k),rm_y(i,jp,k))
             area=scalar_cell_width(gi)*dz
             call record(theta,a*area)
             accepted_integral%y(i,j,k)=theta*a
          enddo
       enddo
    enddo
    do k=1,nz
       kp=modulo(k,nz)+1
       do j=1,zsize(2)
          do i=1,zsize(1)
             gi=i+zstart(1)-1
             if (gi==1) cycle
             a=accepted_integral%z(i,j,k)
             theta=face_theta(a,rp_z(i,j,k),rm_z(i,j,k),rp_z(i,j,kp),rm_z(i,j,kp))
             area=scalar_cell_width(gi)*dy
             call record(theta,a*area)
             accepted_integral%z(i,j,k)=theta*a
          enddo
       enddo
    enddo
    call add_faces(accepted_integral,low_integral(is),1._mytype)
    local(1)=-local(1)
    call MPI_ALLREDUCE(local(1),global(1),1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (ierr/=MPI_SUCCESS) call decomp_2d_abort(1,'Scalar limiter reduction failed')
    call MPI_ALLREDUCE(local(2),global(2),1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    if (ierr/=MPI_SUCCESS) call decomp_2d_abort(1,'Scalar limiter reduction failed')
    call MPI_ALLREDUCE(local(3),global(3),1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (ierr/=MPI_SUCCESS) call decomp_2d_abort(1,'Scalar limiter reduction failed')
    call MPI_ALLREDUCE(count_local,nlimited,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    if (ierr/=MPI_SUCCESS) call decomp_2d_abort(1,'Scalar limiter count failed')
    theta_min=-global(1); removed_l1=global(2); removed_max=global(3)
  contains
    subroutine record(theta,a)
      real(mytype), intent(in) :: theta,a
      real(mytype) :: removed
      if (theta<1._mytype-16._mytype*epsilon(1._mytype)) count_local=count_local+1
      removed=abs((1._mytype-theta)*a)
      local(1)=min(local(1),theta); local(2)=local(2)+removed; local(3)=max(local(3),removed)
    end subroutine record
  end subroutine limit_internal

  ! Only scalar advancement is replaced. int_time calls this after advancing
  ! momentum, but all scalar RHS/flux data were saved from its input velocity.
  subroutine scalar_fct_advance(phi,dphi)
    real(mytype), intent(inout) :: phi(xsize(1),xsize(2),xsize(3),numscalar)
    real(mytype), intent(inout) :: dphi(xsize(1),xsize(2),xsize(3),ntime,numscalar)
    real(mytype) :: bounds(2),identity,global_identity,theta_min,l1,amax
    real(mytype) :: local_mass(3),global_mass(3),volume,residual,tol
    integer :: is,i,j,k,ierr,nlimited
    if (pending_stage/=itr) call decomp_2d_abort(1,'Missing conservative scalar stage RHS')
    do is=1,numscalar
       if (.not.new_scalar(is)) cycle
       if (itr==1) then
          phi(2:nx,:,:,is)=phi(2:nx,:,:,is)+gdt(itr)*dphi(2:nx,:,:,1,is)
       else
          phi(2:nx,:,:,is)=phi(2:nx,:,:,is)+adt(itr)*dphi(2:nx,:,:,1,is)+bdt(itr)*dphi(2:nx,:,:,2,is)
       endif
       dphi(:,:,:,2,is)=dphi(:,:,:,1,is)
       call check_finite(phi(:,:,:,is),'Nonfinite high-order scalar RK update')
       call global_extrema(phi(:,:,:,is),bounds)
       if (nrank==0) write(*,*) 'SCALAR_STAGE step,stage,species,min,max=',itime,itr,is,bounds
       if (itr/=3) cycle
       call scalar_divergence(high_integral(is),div_x)
       identity=maxval(abs(phi(2:nx,:,:,is)-(old(2:nx,:,:,is)-div_x(2:nx,:,:))))
       call MPI_ALLREDUCE(identity,global_identity,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
       if (ierr/=MPI_SUCCESS) call decomp_2d_abort(1,'Scalar RK identity reduction failed')
       if (global_identity>4096._mytype*epsilon(1._mytype)*max(1._mytype,maxval(abs(bounds)))) &
            call decomp_2d_abort(1,'Native scalar RK3 and integrated flux update disagree')
       theta_min=1._mytype; l1=0._mytype; amax=0._mytype; nlimited=0
       if (scalar_scheme(is)=='fct') then
          call limit_internal(is,theta_min,l1,amax,nlimited)
       else
          call faces_copy(high_integral(is),accepted_integral)
       endif
       call scalar_divergence(accepted_integral,div_x)
       phi(2:nx,:,:,is)=old(2:nx,:,:,is)-div_x(2:nx,:,:)
       call check_finite(phi(:,:,:,is),'Nonfinite accepted conservative scalar')
       call global_extrema(phi(:,:,:,is),bounds)
       tol=1024._mytype*epsilon(1._mytype)*max(1._mytype,abs(scalar_fct_lower),abs(scalar_fct_upper))
       if (scalar_scheme(is)=='fct') then
          if (bounds(1)<scalar_fct_lower-tol.or.bounds(2)>scalar_fct_upper+tol) &
               call decomp_2d_abort(1,'FCT accepted scalar outside its independent interval')
       endif
       local_mass=0._mytype
       do k=1,xsize(3)
          do j=1,xsize(2)
             do i=2,nx
                volume=scalar_cell_volume(i)
                local_mass(1)=local_mass(1)+(phi(i,j,k,is)-old(i,j,k,is))*volume
                local_mass(3)=local_mass(3)+abs(old(i,j,k,is))*volume
             enddo
             local_mass(2)=local_mass(2)+(accepted_integral%x(nx,j,k)-accepted_integral%x(1,j,k))*dy*dz
          enddo
       enddo
       call MPI_ALLREDUCE(local_mass,global_mass,3,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
       if (ierr/=MPI_SUCCESS) call decomp_2d_abort(1,'Scalar physical mass reduction failed')
       residual=abs(global_mass(1)+global_mass(2))
       if (residual>8192._mytype*epsilon(1._mytype)*max(1._mytype,global_mass(3))) &
            call decomp_2d_abort(1,'Complete scalar step physical mass budget failed')
       if (nrank==0) then
          write(*,*) 'SCALAR_ACCEPT step,species,scheme=',itime,is,trim(scalar_scheme(is))
          write(*,*) 'SCALAR_ACCEPT min,max,RK_error,mass_error=',bounds,global_identity,residual
          write(*,*) 'SCALAR_MASS delta,physical_outward_integral=',global_mass(1:2)
          write(*,*) 'SCALAR_INTERNAL theta_min,limited_faces,removed_L1,max=',theta_min,nlimited,l1,amax
       endif
       ! Native RK3 ignores history at its first stage. Clear the rejected
       ! high candidate's RHS explicitly, including at full-step checkpoints.
       dphi(:,:,:,:,is)=0._mytype
    enddo
    if (itr==3) pending_stage=0
    sync_scal_needed=.true.
  end subroutine scalar_fct_advance

end module scalar_fct
