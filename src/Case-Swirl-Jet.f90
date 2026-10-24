!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module swirljet

  USE decomp_2d_constants
  USE decomp_2d_mpi
  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  PRIVATE !! All functions/subroutines private by default
  PUBLIC :: init_swirljet, boundary_conditions_swirljet, postprocess_swirljet, &
            visu_swirljet, visu_swirljet_init, scalar_inlet_swirljet

  type :: inlet_table
     integer :: ny = 0, nz = 0
     real(mytype) :: ri = 0, ro = 0
     real(mytype), allocatable :: y(:), z(:), velocity(:,:,:)
  end type inlet_table

  ! Cached local X-pencil inlet plane; loaded once, including after restart.
  real(mytype), allocatable, save :: swirl_velocity(:,:,:)
  logical, parameter :: swirl_derived_output = .false.

contains

  subroutine read_inlet_table(filename, table, status, message)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use, intrinsic :: iso_fortran_env, only : int64
    character(len=*), intent(in) :: filename
    type(inlet_table), intent(out) :: table
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    character(len=:), allocatable :: line
    real(mytype) :: dimensions(2), radii(2)
    integer :: unit, ios, j, k

    status = 1
    message = 'Cannot open swirl inlet file: '//trim(filename)
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) return
    message = 'Invalid swirl inlet header (expected SWIRL_INLET_V1)'
    call read_table_line(unit, line, 32, ios)
    if (ios /= 0) goto 900
    if (trim(line) /= 'SWIRL_INLET_V1') goto 900
    message = 'Invalid swirl inlet dimensions'
    call read_table_row(unit, dimensions, ios)
    if (ios /= 0) goto 900
    if (any(dimensions < 2) .or. any(dimensions > 10000000)) goto 900
    if (any(dimensions /= anint(dimensions))) goto 900
    table%ny = int(dimensions(1))
    table%nz = int(dimensions(2))
    if (table%ny < 2 .or. table%nz < 2) goto 900
    if (int(table%ny,int64)*int(table%nz,int64) > 10000000_int64) goto 900
    message = 'Invalid swirl inlet radii'
    call read_table_row(unit, radii, ios)
    if (ios /= 0) goto 900
    table%ri = radii(1)
    table%ro = radii(2)
    if (.not.ieee_is_finite(table%ri) .or. .not.ieee_is_finite(table%ro)) goto 900
    if (table%ri <= 0 .or. table%ro <= table%ri) goto 900
    allocate(table%y(table%ny), table%z(table%nz), &
         table%velocity(3,table%ny,table%nz), stat=ios)
    message = 'Cannot allocate swirl inlet table'
    if (ios /= 0) goto 900
    message = 'Invalid swirl inlet coordinates'
    call read_table_row(unit, table%y, ios)
    if (ios /= 0) goto 900
    call read_table_row(unit, table%z, ios)
    if (ios /= 0) goto 900
    if (.not.all(ieee_is_finite(table%y)) .or. .not.all(ieee_is_finite(table%z))) goto 900
    if (any(table%y(2:) <= table%y(:table%ny-1))) goto 900
    if (any(table%z(2:) <= table%z(:table%nz-1))) goto 900
    message = 'Invalid or truncated swirl inlet velocity data'
    do k = 1, table%nz
       do j = 1, table%ny
          call read_table_row(unit, table%velocity(:,j,k), ios)
          if (ios /= 0) goto 900
       enddo
    enddo
    if (.not.all(ieee_is_finite(table%velocity))) goto 900
    message = 'Unexpected records after swirl inlet velocity data'
    do
       call read_table_line(unit, line, 256, ios)
       if (ios < 0) exit
       if (ios /= 0) goto 900
       if (len_trim(line) /= 0) goto 900
    enddo
    status = 0
    message = ''
900 close(unit)
  end subroutine read_inlet_table

  ! Read a complete record, rejecting oversized/truncated lines. A negative
  ! status denotes EOF only; other failures must not be mistaken for EOF.
  subroutine read_table_line(unit, line, max_length, status)
    use, intrinsic :: iso_fortran_env, only : iostat_eor, iostat_end
    integer, intent(in) :: unit, max_length
    character(len=:), allocatable, intent(out) :: line
    integer, intent(out) :: status
    character(len=1024) :: chunk
    integer :: ios, received

    line = ''
    status = 1
    do
       read(unit, '(A)', advance='no', size=received, iostat=ios) chunk
       if (received > max_length-len(line)) return
       line = line//chunk(:received)
       if (ios == iostat_eor) then
          status = 0
          return
       endif
       if (ios == iostat_end) then
          status = -1
          return
       endif
       if (ios /= 0) return
    enddo
  end subroutine read_table_line

  ! List-directed input alone accepts null fields, slash termination and repeat
  ! counts. Require exactly one finite numeric token for each output value.
  subroutine read_table_row(unit, values, status)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    integer, intent(in) :: unit
    real(mytype), intent(out) :: values(:)
    integer, intent(out) :: status
    character(len=:), allocatable :: line
    integer :: i, first, last, cursor, length, ios

    values = 0
    call read_table_line(unit, line, 64*size(values), ios)
    status = 1
    if (ios /= 0) return
    length = len_trim(line)
    cursor = 1
    do i = 1, size(values)
       do while (cursor <= length)
          if (line(cursor:cursor) /= ' ' .and. iachar(line(cursor:cursor)) /= 9) exit
          cursor = cursor+1
       enddo
       if (cursor > length) return
       first = cursor
       do while (cursor <= length)
          if (line(cursor:cursor) == ' ' .or. iachar(line(cursor:cursor)) == 9) exit
          cursor = cursor+1
       enddo
       last = cursor-1
       if (verify(line(first:last), '0123456789+-.eEdD') /= 0) return
       read(line(first:last), *, iostat=ios) values(i)
       if (ios /= 0) return
       if (.not.ieee_is_finite(values(i))) return
    enddo
    do while (cursor <= length)
       if (line(cursor:cursor) /= ' ' .and. iachar(line(cursor:cursor)) /= 9) return
       cursor = cursor+1
    enddo
    status = 0
  end subroutine read_table_row

  subroutine sample_inlet_table(table, y, z, velocity, status)
    type(inlet_table), intent(in) :: table
    real(mytype), intent(in) :: y, z
    real(mytype), intent(out) :: velocity(3)
    integer, intent(out) :: status
    integer :: j, k
    real(mytype) :: fy, fz

    velocity = 0
    call bracket(table%y, y, j, fy, status)
    if (status /= 0) return
    call bracket(table%z, z, k, fz, status)
    if (status /= 0) return
    velocity = (1-fy)*(1-fz)*table%velocity(:,j,k) &
         + fy*(1-fz)*table%velocity(:,j+1,k) &
         + (1-fy)*fz*table%velocity(:,j,k+1) &
         + fy*fz*table%velocity(:,j+1,k+1)
  end subroutine sample_inlet_table

  subroutine bracket(axis, value, index, fraction, status)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    real(mytype), intent(in) :: axis(:), value
    integer, intent(out) :: index, status
    real(mytype), intent(out) :: fraction
    integer :: lo, hi, mid
    real(mytype) :: clipped, tolerance

    status = 1
    index = 1
    fraction = 0
    if (.not.ieee_is_finite(value)) return
    tolerance = 64*epsilon(value)*max(1._mytype, abs(axis(1)), abs(axis(size(axis))))
    if (value < axis(1)-tolerance .or. value > axis(size(axis))+tolerance) return
    clipped = max(axis(1), min(value, axis(size(axis))))
    lo = 1
    hi = size(axis)
    do while (hi-lo > 1)
       mid = (lo+hi)/2
       if (clipped < axis(mid)) then
          hi = mid
       else
          lo = mid
       endif
    enddo
    index = lo
    fraction = (clipped-axis(lo))/(axis(lo+1)-axis(lo))
    status = 0
  end subroutine bracket

  ! Validate on every rank before evaluating the annular profile.
  subroutine validate_swirljet_radii()

    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use ibm_param, only : rai, rao
    use mpi, only : MPI_ABORT, MPI_COMM_WORLD

    implicit none

    integer, parameter :: invalid_radii = 1
    integer :: ierror

    if (ieee_is_finite(rai) .and. ieee_is_finite(rao)) then
       if (rai > zero .and. rao > rai) return
    endif

    if (nrank == 0) then
       write(*,*) 'SIMULATION IS STOPPED!'
       write(*,*) 'Swirling jet requires finite radii with 0 < rai < rao in the input file.'
       write(*,*) 'rai = ', rai, ', rao = ', rao
    endif
    call MPI_ABORT(MPI_COMM_WORLD, invalid_radii, ierror)
    stop 1

  end subroutine validate_swirljet_radii

  subroutine load_swirljet_inlet()

    use ibm_param, only : rai, rao
    ! Some Intel MPI bindings provide choice-buffer routines as externals.
    ! Do not require MPI_BCAST/MPI_ALLREDUCE to be module exports via ONLY.
    use mpi

    implicit none

    type(inlet_table) :: table
    integer :: status, global_status, ierr, dims(2), j, k, allocation_status
    character(len=256) :: message
    character(len=1024) :: swirl_inlet_file
    real(mytype) :: y, z, value(3), tolerance

    if (iin < 0 .or. iin > 4 .or. iin == 3) then
       if (nrank == 0) write(*,*) 'Swirling jet supports iin=0/1/2 (analytic) or iin=4 (static table)'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       stop 1
    endif
    if (iin /= 4) return
    if (allocated(swirl_velocity)) return
    if (len_trim(inflowpath) == 0) then
       swirl_inlet_file = 'swirl_inlet.dat'
    else
       swirl_inlet_file = trim(inflowpath)//'/swirl_inlet.dat'
    endif
    status = 0
    if (nrank == 0) then
       if (inflow_noise /= zero) then
          status = 1
          message = 'Static swirl inlet requires inflow_noise=0; use init_noise for initial perturbations'
       else
          call read_inlet_table(trim(swirl_inlet_file), table, status, message)
       endif
       if (status /= 0) write(*,*) trim(message)
    endif
    call MPI_BCAST(status, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (status /= 0) then
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       stop 1
    endif
    if (nrank == 0) dims = [table%ny, table%nz]
    call MPI_BCAST(dims, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    status = 0
    if (nrank /= 0) then
       table%ny = dims(1)
       table%nz = dims(2)
       allocate(table%y(dims(1)), table%z(dims(2)), &
            table%velocity(3,dims(1),dims(2)), stat=allocation_status)
       if (allocation_status /= 0) status = 1
    endif
    call MPI_ALLREDUCE(status, global_status, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    if (global_status /= 0) then
       if (nrank == 0) write(*,*) 'Cannot allocate distributed swirl inlet table'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       stop 1
    endif
    call MPI_BCAST(table%ri, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(table%ro, 1, real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(table%y, dims(1), real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(table%z, dims(2), real_type, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(table%velocity, 3*dims(1)*dims(2), real_type, 0, MPI_COMM_WORLD, ierr)
    tolerance = 64*epsilon(one)*max(one, rao)
    if (abs(table%ri-rai) > tolerance .or. abs(table%ro-rao) > tolerance) then
       if (nrank == 0) write(*,*) 'Swirl inlet table radii do not match rai/rao'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       stop 1
    endif
    allocate(swirl_velocity(3,xsize(2),xsize(3)), stat=allocation_status)
    status = 0
    if (allocation_status /= 0) status = 1
    call MPI_ALLREDUCE(status, global_status, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    if (global_status /= 0) then
       if (nrank == 0) write(*,*) 'Cannot allocate local swirl inlet plane'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       stop 1
    endif
    do k = 1, xsize(3)
       z = (k+xstart(3)-2)*dz-zlz/2
       do j = 1, xsize(2)
          if (istret == 0) then
             y = (j+xstart(2)-2)*dy-yly/2
          else
             y = yp(j+xstart(2)-1)-yly/2
          endif
          call sample_inlet_table(table, y, z, value, allocation_status)
          if (allocation_status /= 0) status = 1
          swirl_velocity(:,j,k) = value
       enddo
    enddo
    call MPI_ALLREDUCE(status, global_status, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    if (global_status /= 0) then
       if (nrank == 0) write(*,*) 'Swirl inlet table does not cover the solver inlet; extrapolation is disabled'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       stop 1
    endif
    if (nrank == 0) write(*,*) 'Static swirl inlet loaded: ', trim(swirl_inlet_file), ' grid=', dims

  end subroutine load_swirljet_inlet

  subroutine boundary_conditions_swirljet (ux,uy,uz,phi)

    USE param
    use ibm_param
    USE variables

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    call inflow (phi)
    call outflow (ux,uy,uz,phi)

    return
  end subroutine boundary_conditions_swirljet

  subroutine inflow (phi)

    USE param
    USE variables
    USE ibm_param

    implicit none

    integer  :: j,k
    real(mytype) :: y, z, r, theta
    real(mytype) :: ratio, um, rm
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    call validate_swirljet_radii()
    call load_swirljet_inlet()

    ratio = rai/rao

    call random_number(bxo)
    call random_number(byo)
    call random_number(bzo)

    do k=1,xsize(3)
       z=(k+xstart(3)-2)*dz-zlz/2.
       do j=1,xsize(2)
          if (istret.eq.0) y=(j+xstart(2)-2)*dy-yly/2.
          if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/2.
          r=sqrt(y**two+z**two)
          rm=r/rao
          theta=atan2(y,z)
          if (iin == 4) then
             bxx1(j,k) = swirl_velocity(1,j,k)
             bxy1(j,k) = swirl_velocity(2,j,k)
             bxz1(j,k) = swirl_velocity(3,j,k)
          elseif (r.gt.rao) then
             bxx1(j,k) = u1/onehundred
             bxy1(j,k) = zero
             bxz1(j,k) = zero
          elseif (r.gt.rai) then
             um=two*(one-rm*rm-log(rm)*(one-ratio*ratio)/log(ratio)) &
               /(one+ratio*ratio+(one-ratio*ratio)/log(ratio))
             bxx1(j,k) = u1*um + u1/onehundred + (bxo(j,k)-half)*inflow_noise
             bxy1(j,k) = u2*rm*cos(theta) + (byo(j,k)-half)*inflow_noise
             bxz1(j,k) = -u2*rm*sin(theta) + (bzo(j,k)-half)*inflow_noise
          else
             bxx1(j,k) = two*(one-r*r/(rai*rai)) + (bxo(j,k)-half)*inflow_noise
             bxy1(j,k) = zero + (byo(j,k)-half)*inflow_noise
             bxz1(j,k) = zero + (bzo(j,k)-half)*inflow_noise
          endif
       enddo
    enddo

    call scalar_inlet_swirljet(phi)

    return
  end subroutine inflow

  subroutine scalar_inlet_swirljet(phi, dphi)

    use ibm_param, only : rai, rao

    real(mytype), intent(inout) :: phi(xsize(1),xsize(2),xsize(3),numscalar)
    real(mytype), intent(inout), optional :: dphi(xsize(1),xsize(2),xsize(3),ntime,numscalar)
    real(mytype), parameter :: delta = 0.1_mytype
    real(mytype) :: y, z, r
    integer :: j, k, is

    if (iscalar.eq.1) then
       do is=1, numscalar
          do k=1,xsize(3)
             z=(k+xstart(3)-2)*dz-zlz/2.
             do j=1,xsize(2)
                if (istret.eq.0) y=(j+xstart(2)-2)*dy-yly/2.
                if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/2.
                r=sqrt(y**two+z**two)
                if (mod(is,2).eq.1) then
                   phi(1,j,k,is)=cp(is)*(one-erf((r-rai)/delta))/two
                else
                   phi(1,j,k,is)=cp(is)*(erf((r-rai)/delta)-erf((r-rao)/delta))/two
                endif
             enddo
          enddo
       enddo

       ! A steady Dirichlet inlet has zero time derivative, including the
       ! stored AB/RK history. Do not retain RHS values from before enforcement.
       if (present(dphi)) dphi(1,:,:,:,:) = zero
    endif

  end subroutine scalar_inlet_swirljet

  subroutine outflow (ux,uy,uz,phi)

    USE param
    USE variables
    USE MPI
    USE ibm_param

    implicit none

    integer :: j,k,code
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx,uxmin,uxmax

    udx=one/dx; udy=one/dy; udz=one/dz; uddx=half/dx; uddy=half/dy; uddz=half/dz

    uxmax=-1609._mytype
    uxmin=1609._mytype
    do k=1,xsize(3)
       do j=1,xsize(2)
          if (ux(nx-1,j,k).gt.uxmax) uxmax=ux(nx-1,j,k)
          if (ux(nx-1,j,k).lt.uxmin) uxmin=ux(nx-1,j,k)
       enddo
    enddo

    call MPI_ALLREDUCE(MPI_IN_PLACE,uxmax,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(MPI_IN_PLACE,uxmin,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)

    if (iin == 4) then
       cx=uxmax*gdt(itr)*udx
    elseif (u1 == zero) then
       cx=(half*(uxmax+uxmin))*gdt(itr)*udx
    elseif (u1 == one) then
       cx=uxmax*gdt(itr)*udx
    elseif (u1 == two) then
       cx=u2*gdt(itr)*udx    !works better
    else
       !cx=(half*(u1+u2))*gdt(itr)*udx
       cx=uxmax*gdt(itr)*udx
    endif

    do k=1,xsize(3)
       do j=1,xsize(2)
          bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
          bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
          bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
       enddo
    enddo

    if (iscalar==1) then
       if (iin == 4) then
          cx=uxmax*gdt(itr)*udx
       elseif (u2==zero) then
          cx=(half*(uxmax+uxmin))*gdt(itr)*udx
       elseif (u2==one) then
          cx=uxmax*gdt(itr)*udx
       elseif (u2==two) then
          cx=u2*gdt(itr)*udx    !works better
       else
          !stop
          cx=uxmax*gdt(itr)*udx
       endif

       do k=1,xsize(3)
          do j=1,xsize(2)
             phi(nx,j,k,:)=phi(nx,j,k,:)-cx*(phi(nx,j,k,:)-phi(nx-1,j,k,:))
          enddo
       enddo
    endif

    if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) &
       write(*,*) "Outflow velocity ux nx=n min max=",real(uxmin,4),real(uxmax,4)

    return
  end subroutine outflow

  subroutine init_swirljet (ux1,uy1,uz1,phi1)

    use navier, only : pre_correc, calc_divu_constraint, solve_poisson, cor_vel
    use var, only : rho1, drho1, ep1, pp3, divu3, px1, py1, pz1
    USE decomp_2d_io
    USE variables
    USE param
    USE ibm_param
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    integer :: i, j, k, ii, is, code
    real(mytype) :: y, z, r, theta
    real(mytype) :: ratio, um, rm, delta

    if (iin.eq.4.and.irestart.eq.0) then
       if (ilmn.or.ivarcoeff.or.iibm.ne.0.or.nclx1.ne.2.or.nclxn.ne.2.or. &
            ncly1.ne.0.or.nclyn.ne.0.or.nclz1.ne.0.or.nclzn.ne.0) then
          call decomp_2d_abort(1, 'Swirl file-inlet initialization requires incompressible flow, '// &
               'no IBM, x Dirichlet and periodic y/z')
       endif
    endif

    call validate_swirljet_radii()
    call load_swirljet_inlet()

    ux1=zero; uy1=zero; uz1=zero

    ratio = rai/rao
    delta = 0.1_mytype

    if (iin.ne.0) then
       call system_clock(count=code)
       if (iin.eq.2) code=0
       call random_seed(size = ii)
       call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))

       call random_number(ux1)
       call random_number(uy1)
       call random_number(uz1)

       !modulation of the random noise
       do k=1,xsize(3)
          z=(k+xstart(3)-2)*dz-zlz/2.
          do j=1,xsize(2)
             if (istret.eq.0) y=(j+xstart(2)-2)*dy-yly/2.
             if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/2.
             r=sqrt(y**two+z**two)
             um=exp(-fifteen*r*r)
             do i=1,xsize(1)
                ux1(i,j,k)=um*init_noise*(ux1(i,j,k)-half)
                uy1(i,j,k)=um*init_noise*(uy1(i,j,k)-half)
                uz1(i,j,k)=um*init_noise*(uz1(i,j,k)-half)
             enddo
          enddo
       enddo
    endif

    !INIT FOR G AND U=MEAN FLOW + NOISE
    do k=1,xsize(3)
       z=(k+xstart(3)-2)*dz-zlz/2.
       do j=1,xsize(2)
          if (istret.eq.0) y=(j+xstart(2)-2)*dy-yly/2.
          if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/2.
          r=sqrt(y**two+z**two)
          rm=r/rao
          theta=atan2(y,z)
          do i=1,xsize(1)
             if (iin == 4) then
                ux1(i,j,k) = ux1(i,j,k) + swirl_velocity(1,j,k)
                uy1(i,j,k) = uy1(i,j,k) + swirl_velocity(2,j,k)
                uz1(i,j,k) = uz1(i,j,k) + swirl_velocity(3,j,k)
             elseif (r.gt.rao) then
                ux1(i,j,k) = ux1(i,j,k) + u1/onehundred
                uy1(i,j,k) = uy1(i,j,k)
                uz1(i,j,k) = uz1(i,j,k)
             elseif (r.gt.rai) then
                um=two*(one-rm*rm-log(rm)*(one-ratio*ratio)/log(ratio)) &
                  /(one+ratio*ratio+(one-ratio*ratio)/log(ratio))
                ux1(i,j,k) = ux1(i,j,k) + u1*um + u1/onehundred
                uy1(i,j,k) = uy1(i,j,k) + u2*rm*cos(theta)
                uz1(i,j,k) = uz1(i,j,k) - u2*rm*sin(theta)
             else
                ux1(i,j,k) = ux1(i,j,k) + two*(one-r*r/(rai*rai))
                uy1(i,j,k) = uy1(i,j,k)
                uz1(i,j,k) = uz1(i,j,k)
             endif
          enddo
       enddo
    enddo

    if (iscalar==1) then
       do k=1,xsize(3)
          z=(k+xstart(3)-2)*dz-zlz/2.
          do j=1,xsize(2)
             if (istret.eq.0) y=(j+xstart(2)-2)*dy-yly/2.
             if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/2.
             r=sqrt(y**two+z**two)
             do i=1,xsize(1)
                do is=1,numscalar
                   if (mod(is,2).eq.1) then
                      phi1(i,j,k,is)=cp(is)*(one-erf((r-rai)/delta))/two
                   else
                      phi1(i,j,k,is)=cp(is)*(erf((r-rai)/delta)-erf((r-rao)/delta))/two
                   endif
                enddo
             enddo
          enddo
       enddo
    endif

    if (iin.eq.4.and.irestart.eq.0) then
       ! Complete u^0 before its first use by either transport equation.
       ! No physical time is advanced and the convective outflow is not called.
       itime = 0; itr = 1
       call inflow(phi1)
       bxxn = bxx1 ! Compatible normal flux for the extruded initial plane.
       bxyn = uy1(nx,:,:); bxzn = uz1(nx,:,:)
       ! Reuse the normal boundary compensation and pressure projection.
       call pre_correc(ux1,uy1,uz1,ep1)
       call calc_divu_constraint(divu3,rho1,phi1)
       call solve_poisson(pp3,px1,py1,pz1,rho1,ux1,uy1,uz1,ep1,drho1,divu3)
       call cor_vel(ux1,uy1,uz1,px1,py1,pz1)
       ! Retain projected u^0, but discard the artificial startup pressure history.
       pp3 = zero; px1 = zero; py1 = zero; pz1 = zero
       dpdyx1 = zero; dpdzx1 = zero
       dpdyxn = zero; dpdzxn = zero
       call postprocess_swirljet(ux1,uy1,uz1,pp3,phi1,ep1)
    endif

#ifdef DEBG
    if (nrank .eq. 0) write(*,*) '# init end ok'
#endif

    return
  end subroutine init_swirljet

  subroutine postprocess_swirljet (ux1,uy1,uz1,pp3,phi1,ep1)

    use mpi
    use var, ONLY : nzmsize
    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

    real(mytype), parameter :: tolerance = 1.e-3_mytype
    real(mytype) :: local(5), global(5), extrema(2), bounds(2), weight, sample_time, value
    integer :: i, j, k, is, code, unit, status, close_status

    if (ipost.eq.0.or.iscalar.eq.0.or.iprocessing.le.0) return
    if (mod(itime,iprocessing).ne.0) return
    if (itime.eq.0) then
       sample_time = t0+real(itime0,mytype)*dt
    else
       sample_time = t
    endif
    ! Separate files keep the new schema distinct from historical statistics.
    call open_scalar_csv_swirljet(unit,'swirl_scalar_bounds.csv', &
         'step,time,scalar,min,max,undershoot_fraction,overshoot_fraction,undershoot_integral,overshoot_integral')
    status = 0
    do is=1,numscalar
       local = zero
       extrema = huge(one)
       do k=1,xsize(3)
          do j=1,xsize(2)
             do i=1,xsize(1)
                weight = one
                if (i.eq.1.or.i.eq.nx) weight = half
                value = phi1(i,j,k,is)
                extrema(1) = min(extrema(1),value)
                extrema(2) = min(extrema(2),-value)
                if (value.lt.-tolerance) local(1) = local(1)+weight
                if (value.gt.one+tolerance) local(2) = local(2)+weight
                local(3) = local(3)+weight*max(-value,zero)
                local(4) = local(4)+weight*max(value-one,zero)
             enddo
          enddo
       enddo
       local(5) = real(nx-1,mytype)*real(xsize(2),mytype)*real(xsize(3),mytype)
       call MPI_ALLREDUCE(local,global,5,real_type,MPI_SUM,MPI_COMM_WORLD,code)
       call MPI_ALLREDUCE(extrema,bounds,2,real_type,MPI_MIN,MPI_COMM_WORLD,code)
       if (nrank.eq.0.and.status.eq.0) &
            write(unit,'(i0,",",es24.16e3,",",i0,6(",",es24.16e3))',iostat=status) &
            itime,sample_time,is,bounds(1),-bounds(2),global(1:2)/global(5),global(3:4)*dx*dy*dz
    enddo
    if (nrank.eq.0) then
       close(unit,iostat=close_status)
       if (status.eq.0) status = close_status
    endif
    call MPI_BCAST(status,1,MPI_INTEGER,0,MPI_COMM_WORLD,code)
    if (status.ne.0) call decomp_2d_abort(1,'Cannot write/close swirl_scalar_bounds.csv')

    ! Bounds every iprocessing steps; costly flux derivatives only every tenth sample.
    if (mod(itime/iprocessing,10).eq.0) &
         call write_scalar_divergence_swirljet(ux1,uy1,uz1,sample_time)

  end subroutine postprocess_swirljet

  subroutine open_scalar_csv_swirljet(unit,filename,header)

    use mpi
    use, intrinsic :: iso_fortran_env, only : int64
    implicit none
    integer, intent(out) :: unit
    character(len=*), intent(in) :: filename, header
    integer :: status, code
    integer(int64) :: bytes

    unit = 0
    status = 0
    if (nrank.eq.0) then
       open(newunit=unit,file=filename,status='unknown',position='append',action='write',iostat=status)
       if (status.eq.0) then
          inquire(unit=unit,size=bytes,iostat=status)
          if (status.eq.0) then
             if (bytes.eq.0) write(unit,'(a)',iostat=status) header
          endif
       endif
    endif
    call MPI_BCAST(status,1,MPI_INTEGER,0,MPI_COMM_WORLD,code)
    if (status.ne.0) call decomp_2d_abort(1,'Cannot open '//filename)

  end subroutine open_scalar_csv_swirljet

  subroutine write_scalar_divergence_swirljet(ux,uy,uz,sample_time)

    use mpi
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    real(mytype), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:), sample_time
    real(mytype), allocatable :: div(:,:,:)
    real(mytype) :: local(2,2), global(2,2), extrema(2,2), bounds(2,2), rms(2), weight
    logical :: last_even
    integer :: i, j, k, is, region, unit, status, close_status, code
    character(len=8), parameter :: regions(2) = ['all     ','interior']

    allocate(div(xsize(1),xsize(2),xsize(3)))
    call open_scalar_csv_swirljet(unit,'swirl_scalar_divergence.csv', &
         'step,time,scalar,region,min,max,rms,volume')
    status = 0
    last_even = .false.
    do is=1,numscalar
       if (is.eq.1.or.(sc_even(is).neqv.last_even)) then
          call scalar_divergence_swirljet(ux,uy,uz,sc_even(is),div)
          local = zero
          extrema = huge(one)
          do k=1,xsize(3)
             do j=1,xsize(2)
                do i=1,xsize(1)
                   weight = one
                   if (i.eq.1.or.i.eq.nx) weight = half
                   local(1,1) = local(1,1)+weight*div(i,j,k)**2
                   local(2,1) = local(2,1)+weight
                   extrema(1,1) = min(extrema(1,1),div(i,j,k))
                   extrema(2,1) = min(extrema(2,1),-div(i,j,k))
                   ! Omit four x planes per end; y/z remain fully periodic.
                   if (i.ge.5.and.i.le.nx-4) then
                      local(1,2) = local(1,2)+div(i,j,k)**2
                      local(2,2) = local(2,2)+one
                      extrema(1,2) = min(extrema(1,2),div(i,j,k))
                      extrema(2,2) = min(extrema(2,2),-div(i,j,k))
                   endif
                enddo
             enddo
          enddo
          call MPI_ALLREDUCE(local,global,4,real_type,MPI_SUM,MPI_COMM_WORLD,code)
          call MPI_ALLREDUCE(extrema,bounds,4,real_type,MPI_MIN,MPI_COMM_WORLD,code)
          do region=1,2
             rms(region) = ieee_value(zero,ieee_quiet_nan)
             if (global(2,region).gt.zero) then
                rms(region) = sqrt(global(1,region)/global(2,region))
             else
                bounds(:,region) = rms(region)
             endif
          enddo
       endif
       last_even = sc_even(is)
       if (nrank.eq.0.and.status.eq.0) then
          do region=1,2
             write(unit,'(i0,",",es24.16e3,",",i0,",",a,4(",",es24.16e3))',iostat=status) &
                  itime,sample_time,is,trim(regions(region)),bounds(1,region),-bounds(2,region), &
                  rms(region),global(2,region)*dx*dy*dz
             if (status.ne.0) exit
          enddo
       endif
    enddo
    if (nrank.eq.0) then
       close(unit,iostat=close_status)
       if (status.eq.0) status = close_status
    endif
    call MPI_BCAST(status,1,MPI_INTEGER,0,MPI_COMM_WORLD,code)
    if (status.ne.0) call decomp_2d_abort(1,'Cannot write/close swirl_scalar_divergence.csv')

  end subroutine write_scalar_divergence_swirljet

  subroutine scalar_divergence_swirljet(ux,uy,uz,even,div)

    use var, only : di1, tb1, ta2, tb2, di2, ta3, tb3, di3
    implicit none
    real(mytype), intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(mytype), intent(out) :: div(:,:,:)
    logical, intent(in) :: even
    ! Match the derivative of u*phi in scalar_transport_eq, including parity
    ! and coefficient arrays. This is D_s.u, not the staggered pressure D_p.u.
    if (even) then
       call derxS(div,ux,di1,sx,ffxS,fsxS,fwxS,xsize(1),xsize(2),xsize(3),0,zero)
    else
       call derxS(div,ux,di1,sx,ffxpS,fsxpS,fwxpS,xsize(1),xsize(2),xsize(3),1,zero)
    endif
    call transpose_x_to_y(uy,ta2)
    if (even) then
       call deryS(tb2,ta2,di2,sy,ffyS,fsyS,fwyS,ppy,ysize(1),ysize(2),ysize(3),0,zero)
    else
       call deryS(tb2,ta2,di2,sy,ffypS,fsypS,fwypS,ppy,ysize(1),ysize(2),ysize(3),1,zero)
    endif
    call transpose_y_to_x(tb2,tb1)
    div = div+tb1
    call transpose_x_to_y(uz,ta2)
    call transpose_y_to_z(ta2,ta3)
    if (even) then
       call derzS(tb3,ta3,di3,sz,ffzS,fszS,fwzS,zsize(1),zsize(2),zsize(3),0,zero)
    else
       call derzS(tb3,ta3,di3,sz,ffzpS,fszpS,fwzpS,zsize(1),zsize(2),zsize(3),1,zero)
    endif
    call transpose_z_to_y(tb3,tb2)
    call transpose_y_to_x(tb2,tb1)
    div = div+tb1

  end subroutine scalar_divergence_swirljet

  subroutine visu_swirljet_init (visu_initialised)

    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D

    implicit none

    logical, intent(out) :: visu_initialised

    if (.not.swirl_derived_output) then
       visu_initialised = .true.
       return
    endif

    call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "critq", 1, 0, output2D, mytype)

    visu_initialised = .true.

  end subroutine visu_swirljet_init

  subroutine visu_swirljet (ux1, uy1, uz1, pp3, phi1, ep1, num)

    use var, only : ux2, uy2, uz2, ux3, uy3, uz3
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2
    USE var, only : ta3,tb3,tc3,td3,te3,tf3,di3
    use var, ONLY : nxmsize, nymsize, nzmsize
    use visu, only : write_field
    use ibm_param, only : ubcx,ubcy,ubcz

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    integer, intent(in) :: num

    integer :: i, j, k
    real(mytype) :: y, z, theta, tmp
    real(mytype) :: eig1, eig2, eig3
    real(mytype) :: m11, m12, m13, m22, m23, m33, trace, det, b, q, r
    real(mytype), allocatable :: ur(:,:,:), utheta(:,:,:)

    if (.not.swirl_derived_output) return
    allocate(ur(xsize(1),xsize(2),xsize(3)),utheta(xsize(1),xsize(2),xsize(3)))

    ! Write vorticity as an example of post processing

    ! Perform communications if needed
    if (sync_vel_needed) then
      call transpose_x_to_y(ux1,ux2)
      call transpose_x_to_y(uy1,uy2)
      call transpose_x_to_y(uz1,uz2)
      call transpose_y_to_z(ux2,ux3)
      call transpose_y_to_z(uy2,uy3)
      call transpose_y_to_z(uz2,uz3)
      sync_vel_needed = .false.
    endif

    !x-derivatives
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    !y-derivatives
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    !!z-derivatives
    call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
    call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
    call derz (tc3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)
    !!all back to x-pencils
    call transpose_z_to_y(ta3,td2)
    call transpose_z_to_y(tb3,te2)
    call transpose_z_to_y(tc3,tf2)
    call transpose_y_to_x(td2,tg1)
    call transpose_y_to_x(te2,th1)
    call transpose_y_to_x(tf2,ti1)
    call transpose_y_to_x(ta2,td1)
    call transpose_y_to_x(tb2,te1)
    call transpose_y_to_x(tc2,tf1)
    !du/dx=ta1 du/dy=td1 and du/dz=tg1
    !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
    !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1

    !VORTICITY FIELD
    di1 = zero
    di1(:,:,:)=sqrt(  (tf1(:,:,:)-th1(:,:,:))**2 &
                    + (tg1(:,:,:)-tc1(:,:,:))**2 &
                    + (tb1(:,:,:)-td1(:,:,:))**2)
    call write_field(di1, ".", "vort", num, flush = .true.) ! Reusing temporary array, force flush

    !Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
    di1 = zero
    di1(:,:,:) = - half*(ta1(:,:,:)**2 + te1(:,:,:)**2 + ti1(:,:,:)**2) &
                 - td1(:,:,:) * tb1(:,:,:) &
                 - tg1(:,:,:) * tc1(:,:,:) &
                 - th1(:,:,:) * tf1(:,:,:)
    call write_field(di1, ".", "critq", num, flush = .true.) ! Reusing temporary array, force flush

    do k=1,xsize(3)
       z=(k+xstart(3)-2)*dz-zlz/2.
       do j=1,xsize(2)
          if (istret.eq.0) y=(j+xstart(2)-2)*dy-yly/2.
          if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/2.
          theta=atan2(y,z)
          do i=1,xsize(1)
             ur(i,j,k) = uz1(i,j,k)*cos(theta) + uy1(i,j,k)*sin(theta)
             utheta(i,j,k) = -uz1(i,j,k)*sin(theta) + uy1(i,j,k)*cos(theta)
          enddo
       enddo
    enddo
    call write_field(ur, ".", "ur", num)
    call write_field(utheta, ".", "utheta", num)

  end subroutine visu_swirljet

end module swirljet
