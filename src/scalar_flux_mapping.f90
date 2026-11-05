!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module scalar_flux_mapping

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : real_type, decomp_2d_abort, nrank
  use decomp_2d, only : xsize, ysize, zsize, xstart, ystart, zstart, &
       transpose_x_to_y, transpose_y_to_z
  use iso_fortran_env, only : int64, output_unit
  use variables, only : nx, ny, nz
  use mpi
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use param, only : dx, dy, dz, itr, itime, scalar_mapping_tol, &
       scalar_mapping_flux_tol, scalar_mapping_maxiter, scalar_mapping_max_correction, &
       scalar_mapping_reference_velocity, scalar_mapping_local_velocity_floor, scalar_mapping_near_x
  use scalar_flux, only : scalar_faces, faces_allocate, scalar_raw_velocity, &
       scalar_divergence, scalar_cell_volume, scalar_cell_width

  implicit none
  private

  public :: scalar_map_velocity
  public :: scalar_map_last_iterations, scalar_map_last_max_du
  public :: scalar_map_last_mean_du, scalar_map_last_mean_abs_du
  public :: scalar_map_last_mean_du_xyz
  public :: scalar_map_last_divergence, scalar_map_last_raw_divergence
  public :: scalar_map_last_compatibility, scalar_map_last_seconds

  integer, save :: scalar_map_last_iterations = 0
  real(mytype), save :: scalar_map_last_max_du = 0._mytype
  real(mytype), save :: scalar_map_last_mean_du = 0._mytype
  real(mytype), save :: scalar_map_last_mean_abs_du = 0._mytype
  real(mytype), save :: scalar_map_last_mean_du_xyz(3) = 0._mytype
  real(mytype), save :: scalar_map_last_divergence = 0._mytype
  real(mytype), save :: scalar_map_last_raw_divergence = 0._mytype
  real(mytype), save :: scalar_map_last_compatibility = 0._mytype
  real(mytype), save :: scalar_map_last_seconds = 0._mytype

  ! All vectors have one local pencil partition, including the inactive inlet
  ! node. No rank stores the whole physical field or reduces a field globally.
  logical, save :: initialized = .false.
  type(scalar_faces), save :: gradient
  real(mytype), allocatable, save :: lambda_x(:,:,:), residual_x(:,:,:)
  real(mytype), allocatable, save :: direction_x(:,:,:), image_x(:,:,:)
  real(mytype), allocatable, save :: preconditioned_x(:,:,:)
  real(mytype), allocatable, save :: lambda_y(:,:,:), lambda_z(:,:,:)
  real(mytype), allocatable, save :: volume(:), diagonal(:)
  real(mytype), save :: global_volume

contains

  subroutine mapping_initialize()

    integer :: i, ierr
    real(mytype) :: local_volume, x_diagonal

    if (initialized) return
    call faces_allocate(gradient)
    allocate(lambda_x(xsize(1),xsize(2),xsize(3)))
    allocate(residual_x(xsize(1),xsize(2),xsize(3)))
    allocate(direction_x(xsize(1),xsize(2),xsize(3)))
    allocate(image_x(xsize(1),xsize(2),xsize(3)))
    allocate(preconditioned_x(xsize(1),xsize(2),xsize(3)))
    allocate(lambda_y(ysize(1),ysize(2),ysize(3)))
    allocate(lambda_z(zsize(1),zsize(2),zsize(3)))
    allocate(volume(xsize(1)),diagonal(xsize(1)))
    volume(1) = 0._mytype
    diagonal(1) = 1._mytype
    do i = 2, xsize(1)
       volume(i) = scalar_cell_volume(i)
       x_diagonal = 0._mytype
       if (i > 2) x_diagonal = x_diagonal + 1._mytype/dx
       if (i < xsize(1)) x_diagonal = x_diagonal + 1._mytype/dx
       diagonal(i) = x_diagonal/scalar_cell_width(i) + 2._mytype/dy**2 + 2._mytype/dz**2
    enddo
    ! These are diagonals of div(G), not D2_iSVV. The first active volume
    ! covers the prescribed inlet gap; the outlet volume has half width.
    local_volume = sum(volume)*real(xsize(2)*xsize(3),mytype)
    call MPI_ALLREDUCE(local_volume,global_volume,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)
    if (global_volume <= 0._mytype) call mapping_abort('Nonpositive active scalar volume')
    initialized = .true.

  end subroutine mapping_initialize

  ! Geometric face-gradient space. It has a constant nullspace only; unlike
  ! compact nodal transverse derivatives it has no even-grid Nyquist kernel.
  ! Positive face density leaves the left cell. Physical x faces stay zero.
  subroutine face_gradient(field)

    real(mytype), intent(in) :: field(xsize(1),xsize(2),xsize(3))
    integer :: i, j, k, jp, kp

    gradient%x = 0._mytype
    do k = 1, xsize(3)
       do j = 1, xsize(2)
          do i = 2, xsize(1)-1
             gradient%x(i,j,k) = (field(i,j,k)-field(i+1,j,k))/dx
          enddo
       enddo
    enddo
    call transpose_x_to_y(field,lambda_y)
    gradient%y = 0._mytype
    do k = 1, ysize(3)
       do j = 1, ysize(2)
          jp = modulo(j,ysize(2))+1
          do i = 1, ysize(1)
             if (ystart(1)+i-1 == 1) cycle
             gradient%y(i,j,k) = (lambda_y(i,j,k)-lambda_y(i,jp,k))/dy
          enddo
       enddo
    enddo
    call transpose_y_to_z(lambda_y,lambda_z)
    gradient%z = 0._mytype
    do k = 1, zsize(3)
       kp = modulo(k,zsize(3))+1
       do j = 1, zsize(2)
          do i = 1, zsize(1)
             if (zstart(1)+i-1 == 1) cycle
             gradient%z(i,j,k) = (lambda_z(i,j,k)-lambda_z(i,j,kp))/dz
          enddo
       enddo
    enddo

  end subroutine face_gradient

  function weighted_dot(a,b) result(value)

    real(mytype), intent(in) :: a(xsize(1),xsize(2),xsize(3))
    real(mytype), intent(in) :: b(xsize(1),xsize(2),xsize(3))
    real(mytype) :: value, local_value
    integer :: i, ierr

    local_value = 0._mytype
    do i = 2, xsize(1)
       local_value = local_value + volume(i)*sum(a(i,:,:)*b(i,:,:))
    enddo
    call MPI_ALLREDUCE(local_value,value,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)

  end function weighted_dot

  function weighted_mean(a) result(value)

    real(mytype), intent(in) :: a(xsize(1),xsize(2),xsize(3))
    real(mytype) :: value, local_value
    integer :: i, ierr

    local_value = 0._mytype
    do i = 2, xsize(1)
       local_value = local_value + volume(i)*sum(a(i,:,:))
    enddo
    call MPI_ALLREDUCE(local_value,value,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)
    value = value/global_volume

  end function weighted_mean

  function global_max_abs(a) result(value)

    real(mytype), intent(in) :: a(xsize(1),xsize(2),xsize(3))
    real(mytype) :: value, local_value
    integer :: ierr

    local_value = maxval(abs(a(2:xsize(1),:,:)))
    call MPI_ALLREDUCE(local_value,value,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)

  end function global_max_abs

  subroutine zero_mean(a)

    real(mytype), intent(inout) :: a(xsize(1),xsize(2),xsize(3))
    real(mytype) :: average

    average = weighted_mean(a)
    a(2:xsize(1),:,:) = a(2:xsize(1),:,:) - average
    a(1,:,:) = 0._mytype

  end subroutine zero_mean

  subroutine precondition()

    integer :: i

    preconditioned_x(1,:,:) = 0._mytype
    do i = 2, xsize(1)
       preconditioned_x(i,:,:) = residual_x(i,:,:)/diagonal(i)
    enddo
    ! Project the Jacobi inverse into the compatible zero-mean subspace.
    call zero_mean(preconditioned_x)

  end subroutine precondition

  subroutine correction_diagnostics(ux,uy,uz,raw)

    real(mytype), intent(in) :: ux(xsize(1),xsize(2),xsize(3))
    real(mytype), intent(in) :: uy(xsize(1),xsize(2),xsize(3))
    real(mytype), intent(in) :: uz(xsize(1),xsize(2),xsize(3))
    type(scalar_faces), intent(in) :: raw
    real(mytype) :: local(10), global(10), weight, speed, scale
    real(mytype) :: region_sum(4,2), global_region_sum(4,2)
    real(mytype) :: region_max(3,2), global_region_max(3,2)
    real(mytype) :: mean_abs(2), rms(2), floored_fraction(2), floored_weight_fraction(2)
    real(mytype) :: local_best_real(9), selected_real(9), point(3)
    integer(int64) :: region_count(2,2), global_region_count(2,2), best_key, candidate_key, selected_key
    integer :: i, j, k, jp, kp, gi, gj, gk, ierr, region, invalid, global_invalid
    integer :: local_best_integer(6), selected_integer(6), owner_candidate, owner

    ! The potential transposes are no longer needed after gradient is formed.
    ! Reuse their storage for raw nodal full-speed magnitudes, not just the
    ! normal component that may vanish on a strongly tangential face.
    do k = 1, xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)
             scale = max(abs(ux(i,j,k)),abs(uy(i,j,k)),abs(uz(i,j,k)))
             image_x(i,j,k) = 0._mytype
             if (scale > 0._mytype) image_x(i,j,k) = scale*sqrt((ux(i,j,k)/scale)**2 &
                  +(uy(i,j,k)/scale)**2+(uz(i,j,k)/scale)**2)
          enddo
       enddo
    enddo
    invalid = 0
    if (.not.all(ieee_is_finite(image_x)) .or. .not.all(ieee_is_finite(gradient%x)) .or. &
         .not.all(ieee_is_finite(gradient%y)) .or. .not.all(ieee_is_finite(gradient%z))) invalid = 1
    call MPI_ALLREDUCE(invalid,global_invalid,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)
    if (global_invalid /= 0) call mapping_abort('Nonfinite mapping correction or diagnostic velocity scale')
    call transpose_x_to_y(image_x,lambda_y)
    call transpose_y_to_z(lambda_y,lambda_z)

    local = 0._mytype
    region_sum = 0._mytype
    region_max = 0._mytype
    region_count = 0_int64
    best_key = huge(best_key)
    local_best_integer = 0
    local_best_real = 0._mytype
    do k = 1, xsize(3)
       gk = xstart(3)+k-1
       do j = 1, xsize(2)
          gj = xstart(2)+j-1
          do i = 2, xsize(1)-1
             weight = dx*dy*dz
             speed = 0.5_mytype*(image_x(i,j,k)+image_x(i+1,j,k))
             point = [(real(i,mytype)-0.5_mytype)*dx,real(gj-1,mytype)*dy,real(gk-1,mytype)*dz]
             call record(-gradient%x(i,j,k),raw%x(i,j,k),speed,weight,1,i,gj,gk,point)
          enddo
       enddo
    enddo
    do k = 1, ysize(3)
       gk = ystart(3)+k-1
       do j = 1, ysize(2)
          jp = modulo(j,ysize(2))+1
          do i = 1, ysize(1)
             gi = ystart(1)+i-1
             if (gi == 1) cycle
             weight = scalar_cell_volume(gi)
             speed = 0.5_mytype*(lambda_y(i,j,k)+lambda_y(i,jp,k))
             point = [real(gi-1,mytype)*dx,(real(j,mytype)-0.5_mytype)*dy,real(gk-1,mytype)*dz]
             call record(-gradient%y(i,j,k),raw%y(i,j,k),speed,weight,2,gi,j,gk,point)
          enddo
       enddo
    enddo
    do k = 1, zsize(3)
       kp = modulo(k,zsize(3))+1
       do j = 1, zsize(2)
          gj = zstart(2)+j-1
          do i = 1, zsize(1)
             gi = zstart(1)+i-1
             if (gi == 1) cycle
             weight = scalar_cell_volume(gi)
             speed = 0.5_mytype*(lambda_z(i,j,k)+lambda_z(i,j,kp))
             point = [real(gi-1,mytype)*dx,real(gj-1,mytype)*dy,(real(k,mytype)-0.5_mytype)*dz]
             call record(-gradient%z(i,j,k),raw%z(i,j,k),speed,weight,3,gi,gj,k,point)
          enddo
       enddo
    enddo
    call MPI_ALLREDUCE(local(1),global(1),1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)
    call MPI_ALLREDUCE(local(2:10),global(2:10),9,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)
    scalar_map_last_max_du = global(1)
    scalar_map_last_mean_du = global(2)/global(4)
    scalar_map_last_mean_abs_du = global(3)/global(4)
    scalar_map_last_mean_du_xyz = global(5:7)/global(8:10)

    call MPI_ALLREDUCE(region_sum,global_region_sum,8,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)
    call MPI_ALLREDUCE(region_max,global_region_max,6,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)
    call MPI_ALLREDUCE(region_count,global_region_count,4,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)
    mean_abs = 0._mytype
    rms = 0._mytype
    floored_fraction = 0._mytype
    floored_weight_fraction = 0._mytype
    do region = 1, 2
       if (global_region_sum(1,region) > 0._mytype) then
          mean_abs(region) = global_region_sum(2,region)/global_region_sum(1,region)
          rms(region) = sqrt(global_region_sum(3,region)/global_region_sum(1,region))
          floored_weight_fraction(region) = global_region_sum(4,region)/global_region_sum(1,region)
       endif
       if (global_region_count(1,region) > 0_int64) floored_fraction(region) = &
            real(global_region_count(2,region),mytype)/real(global_region_count(1,region),mytype)
    enddo

    ! All corrections are owned exactly once in their directional pencils.
    ! Among ties at the global maximum choose the smallest global face key,
    ! then its owner rank. Only a short metadata record is broadcast.
    candidate_key = huge(candidate_key)
    if (local(1) == global(1)) candidate_key = best_key
    call MPI_ALLREDUCE(candidate_key,selected_key,1,MPI_INTEGER8,MPI_MIN,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)
    owner_candidate = huge(owner_candidate)
    if (candidate_key == selected_key) owner_candidate = nrank
    call MPI_ALLREDUCE(owner_candidate,owner,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)
    selected_integer = 0
    selected_real = 0._mytype
    if (nrank == owner) then
       selected_integer = local_best_integer
       selected_real = local_best_real
    endif
    call MPI_BCAST(selected_integer,6,MPI_INTEGER,owner,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)
    call MPI_BCAST(selected_real,9,real_type,owner,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)
    if (nrank == 0) then
       write(*,*) 'SCALAR_MAP correction_diagnostics step,stage=',itime,itr
       write(*,*) 'SCALAR_MAP reference_velocity,local_velocity_floor,near_x=', &
            scalar_mapping_reference_velocity,scalar_mapping_local_velocity_floor,scalar_mapping_near_x
       write(*,*) 'SCALAR_MAP max_location_ijk,component,rank,floored=',selected_integer
       write(*,*) 'SCALAR_MAP max_location_xyz_solver_origin=',selected_real(1:3)
       write(*,*) 'SCALAR_MAP max_signed_delta_u,raw_normal_velocity,raw_local_speed,used_local_scale=',selected_real(4:7)
       write(*,*) 'SCALAR_MAP max_location_delta_over_reference,delta_over_local=',selected_real(8:9)
       write(*,*) 'SCALAR_MAP max_abs_delta_over_reference,max_local_relative=', &
            global(1)/scalar_mapping_reference_velocity,maxval(global_region_max(3,:))
       do region = 1, 2
          write(*,*) 'SCALAR_MAP region_index(1=near_x_le_cutoff,2=far),face_count,floored_count=', &
               region,global_region_count(:,region)
          write(*,*) 'SCALAR_MAP region_abs_mean,RMS,max,max_over_reference,max_over_local=', &
               mean_abs(region),rms(region),global_region_max(:,region)
          write(*,*) 'SCALAR_MAP region_floored_fraction,weighted_floored_fraction=', &
               floored_fraction(region),floored_weight_fraction(region)
       enddo
    endif

  contains

    subroutine record(delta,raw_normal,raw_speed,face_weight,component,global_i,global_j,global_k,xyz)

      real(mytype), intent(in) :: delta,raw_normal,raw_speed,face_weight,xyz(3)
      integer, intent(in) :: component,global_i,global_j,global_k
      real(mytype) :: magnitude,denominator,relative_reference,relative_local
      integer(int64) :: key
      integer :: r,floored

      magnitude = abs(delta)
      denominator = max(raw_speed,scalar_mapping_local_velocity_floor)
      relative_reference = magnitude/scalar_mapping_reference_velocity
      relative_local = magnitude/denominator
      floored = 0
      if (raw_speed < scalar_mapping_local_velocity_floor) floored = 1
      r = 2
      if (xyz(1) <= scalar_mapping_near_x) r = 1
      region_count(1,r) = region_count(1,r)+1_int64
      region_count(2,r) = region_count(2,r)+int(floored,int64)
      region_sum(1,r) = region_sum(1,r)+face_weight
      region_sum(2,r) = region_sum(2,r)+face_weight*magnitude
      region_sum(3,r) = region_sum(3,r)+face_weight*delta**2
      if (floored /= 0) region_sum(4,r) = region_sum(4,r)+face_weight
      region_max(1,r) = max(region_max(1,r),magnitude)
      region_max(2,r) = max(region_max(2,r),relative_reference)
      region_max(3,r) = max(region_max(3,r),relative_local)
      local(2) = local(2)+face_weight*delta
      local(3) = local(3)+face_weight*magnitude
      local(4) = local(4)+face_weight
      local(4+component) = local(4+component)+face_weight*delta
      local(7+component) = local(7+component)+face_weight
      key = (((int(component-1,int64)*int(nx,int64)+int(global_i-1,int64))*int(ny,int64) &
           +int(global_j-1,int64))*int(nz,int64))+int(global_k,int64)
      if (magnitude > local(1) .or. (magnitude == local(1).and.key < best_key)) then
         local(1) = magnitude
         best_key = key
         local_best_integer = [global_i,global_j,global_k,component,nrank,floored]
         local_best_real = [xyz,delta,raw_normal,raw_speed,denominator,delta/scalar_mapping_reference_velocity, &
              delta/denominator]
      endif

    end subroutine record

  end subroutine correction_diagnostics

  ! Called once per real RHS stage and shared by all passive species. This
  ! changes only auxiliary face velocities, never the nodal momentum field.
  ! The volume-weighted Laplacian is SPD on the zero-mean active subspace:
  ! <lambda,L lambda>_V = sum_internal (area/distance)*(delta lambda)^2.
  subroutine scalar_map_velocity(ux,uy,uz,mapped)

    real(mytype), intent(in) :: ux(xsize(1),xsize(2),xsize(3))
    real(mytype), intent(in) :: uy(xsize(1),xsize(2),xsize(3))
    real(mytype), intent(in) :: uz(xsize(1),xsize(2),xsize(3))
    type(scalar_faces), intent(inout) :: mapped
    real(mytype) :: rz, rz_new, denominator, alpha, residual_norm
    real(mytype) :: local_boundary_min, boundary_min, start_time, elapsed
    integer :: iteration, invalid, global_invalid, ierr

    call mapping_initialize()
    call faces_allocate(mapped)
    start_time = real(MPI_WTIME(),mytype)
    scalar_map_last_iterations = 0
    scalar_map_last_max_du = 0._mytype
    scalar_map_last_mean_du = 0._mytype
    scalar_map_last_mean_abs_du = 0._mytype
    scalar_map_last_mean_du_xyz = 0._mytype

    invalid = 0
    if (.not.all(ieee_is_finite(ux)) .or. .not.all(ieee_is_finite(uy)) .or. &
         .not.all(ieee_is_finite(uz))) invalid = 1
    call MPI_ALLREDUCE(invalid,global_invalid,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)
    if (global_invalid /= 0) call mapping_abort('Nonfinite real-stage velocity')

    call scalar_raw_velocity(ux,uy,uz,mapped)
    invalid = 0
    if (.not.all(ieee_is_finite(mapped%x)) .or. .not.all(ieee_is_finite(mapped%y)) .or. &
         .not.all(ieee_is_finite(mapped%z))) invalid = 1
    call MPI_ALLREDUCE(invalid,global_invalid,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)
    if (global_invalid /= 0) call mapping_abort('Nonfinite reconstructed scalar mass velocity')
    local_boundary_min = min(minval(mapped%x(1,:,:)),minval(mapped%x(xsize(1),:,:)))
    call MPI_ALLREDUCE(local_boundary_min,boundary_min,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)
    if (boundary_min < -scalar_mapping_tol) then
       if (nrank == 0) write(*,*) 'SCALAR_MAP unsupported physical backflow:',boundary_min
       call mapping_abort('Scalar auxiliary mapping does not support physical backflow')
    endif

    call scalar_divergence(mapped,residual_x)
    scalar_map_last_raw_divergence = global_max_abs(residual_x)
    scalar_map_last_compatibility = weighted_mean(residual_x)
    if (abs(scalar_map_last_compatibility) > scalar_mapping_flux_tol) then
       if (nrank == 0) write(*,*) 'SCALAR_MAP incompatible fixed physical flux: mean div, net flow=', &
            scalar_map_last_compatibility,scalar_map_last_compatibility*global_volume
       call mapping_abort('Fixed physical mass fluxes are incompatible; no automatic outlet repair')
    endif
    ! An allowed net-flow tolerance cannot override the requested local
    ! divergence tolerance. An internal correction cannot remove this mode.
    if (abs(scalar_map_last_compatibility) > scalar_mapping_tol) then
       call mapping_abort('Fixed physical mean divergence exceeds scalar_mapping_tol')
    endif
    lambda_x = 0._mytype
    call zero_mean(residual_x)
    call precondition()
    direction_x = preconditioned_x
    rz = weighted_dot(residual_x,preconditioned_x)
    residual_norm = global_max_abs(residual_x)
    iteration = 0
    do while (residual_norm > 0.5_mytype*scalar_mapping_tol)
       if (iteration >= scalar_mapping_maxiter) then
          if (nrank == 0) write(*,*) 'SCALAR_MAP CG iterations/residual=',iteration,residual_norm
          call mapping_abort('Scalar auxiliary mapping did not converge')
       endif
       call face_gradient(direction_x)
       call scalar_divergence(gradient,image_x)
       denominator = weighted_dot(direction_x,image_x)
       if (.not.ieee_is_finite(denominator) .or. denominator <= tiny(denominator) .or. &
            .not.ieee_is_finite(rz) .or. rz <= 0._mytype) then
          call mapping_abort('Scalar auxiliary mapping CG lost positive definiteness')
       endif
       alpha = rz/denominator
       lambda_x = lambda_x + alpha*direction_x
       residual_x = residual_x - alpha*image_x
       call zero_mean(residual_x)
       iteration = iteration+1
       residual_norm = global_max_abs(residual_x)
       if (residual_norm <= 0.5_mytype*scalar_mapping_tol) exit
       call precondition()
       rz_new = weighted_dot(residual_x,preconditioned_x)
       direction_x = preconditioned_x + (rz_new/rz)*direction_x
       rz = rz_new
    enddo
    scalar_map_last_iterations = iteration
    call face_gradient(lambda_x)
    call correction_diagnostics(ux,uy,uz,mapped)
    if (scalar_map_last_max_du > scalar_mapping_max_correction) then
       if (nrank == 0) then
          write(*,*) 'Scalar auxiliary velocity correction exceeds configured budget'
          write(*,*) 'SCALAR_MAP internal correction/budget=', &
               scalar_map_last_max_du,scalar_mapping_max_correction
       endif
       ! Preserve the root's diagnostic context before any rank aborts.
       flush(output_unit)
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call check_mpi(ierr)
       call mapping_abort('Scalar auxiliary velocity correction exceeds configured budget')
    endif
    mapped%x = mapped%x-gradient%x
    mapped%y = mapped%y-gradient%y
    mapped%z = mapped%z-gradient%z
    call scalar_divergence(mapped,image_x)
    scalar_map_last_divergence = global_max_abs(image_x)
    if (.not.ieee_is_finite(scalar_map_last_divergence) .or. &
         scalar_map_last_divergence > 2._mytype*scalar_mapping_tol) then
       if (nrank == 0) write(*,*) 'SCALAR_MAP final divergence=',scalar_map_last_divergence
       call mapping_abort('Scalar auxiliary face divergence failed final verification')
    endif
    elapsed = real(MPI_WTIME(),mytype)-start_time
    call MPI_ALLREDUCE(elapsed,scalar_map_last_seconds,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    call check_mpi(ierr)
    if (nrank == 0) then
       write(*,*) 'SCALAR_MAP step,stage,iterations=',itime,itr,scalar_map_last_iterations
       write(*,*) 'SCALAR_MAP raw_div,mapped_div,mean_div=',scalar_map_last_raw_divergence, &
            scalar_map_last_divergence,scalar_map_last_compatibility
       write(*,*) 'SCALAR_MAP max_du,mean_du,mean_abs_du,seconds=',scalar_map_last_max_du, &
            scalar_map_last_mean_du,scalar_map_last_mean_abs_du,scalar_map_last_seconds
       write(*,*) 'SCALAR_MAP mean_du_xyz=',scalar_map_last_mean_du_xyz
       write(*,*) 'SCALAR_MAP physical_max_du,physical_net_change=',0._mytype,0._mytype
       write(*,*) 'SCALAR_MAP physical_outlet_delta_modes=',0._mytype,0._mytype,0._mytype,0._mytype
    endif

  end subroutine scalar_map_velocity

  subroutine check_mpi(ierr)

    integer, intent(in) :: ierr

    if (ierr /= MPI_SUCCESS) call mapping_abort('MPI failure in scalar auxiliary mapping')

  end subroutine check_mpi

  subroutine mapping_abort(message)

    character(len=*), intent(in) :: message

    call decomp_2d_abort(__FILE__,__LINE__,1,message)

  end subroutine mapping_abort

end module scalar_flux_mapping
