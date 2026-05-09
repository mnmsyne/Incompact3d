!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module walljet

  USE decomp_2d_constants
  USE decomp_2d_mpi
  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  character(len=*), parameter :: io_jet = "io-jet"

  PRIVATE !! All functions/subroutines private by default
  PUBLIC :: init_walljet, &
            boundary_conditions_walljet, &
            momentum_forcing_walljet, &
            postprocess_walljet, &
            visu_walljet, &
            visu_walljet_init

contains

   subroutine boundary_conditions_walljet(ux, uy, uz, phi)

      use param
      use variables

      implicit none

      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux, uy, uz
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi
      integer :: i, k
      real(mytype) :: ramp, top_v_eff

      ramp = one
      if (walljet_ramp_tau > tiny(one)) then
         ramp = one - exp(-max(zero, real(itime - 1, mytype) * dt) / walljet_ramp_tau)
      end if
      top_v_eff = one - ramp * (one - max(zero, min(one, walljet_top_v_relax)))

      call inflow(phi)
      call outflow(ux, uy, uz, phi)

      ! Bottom wall (no-slip)
      if (ncly1 == 2) then
         do k = 1, xsize(3)
            do i = 1, xsize(1)
               byx1(i, k) = zero
               byy1(i, k) = zero
               byz1(i, k) = zero
            end do
         end do
      end if

      ! Top boundary (Dirichlet uses interior extrapolation)
      if (nclyn == 2) then
         do k = 1, xsize(3)
            do i = 1, xsize(1)
               byxn(i, k) = ux(i, xsize(2) - 1, k)
               byyn(i, k) = top_v_eff * uy(i, xsize(2) - 1, k)
               byzn(i, k) = uz(i, xsize(2) - 1, k)
            end do
         end do
      end if

      return
   end subroutine boundary_conditions_walljet
  !############################################################################
  
   subroutine inflow(phi)

      use param
      use variables
      use ibm_param

      implicit none

      integer :: j, k, is
      real(mytype) :: y, um1, um2, um, ramp, ubase
      real(mytype), dimension(xsize(2), xsize(3)) :: bcx, bcy, bcz
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi

      ramp = one
      if (walljet_ramp_tau > tiny(one)) then
         ramp = one - exp(-max(zero, real(itime - 1, mytype) * dt) / walljet_ramp_tau)
      end if

      call random_number(bcx)
      call random_number(bcy)
      call random_number(bcz)

      do k = 1, xsize(3)
         do j = 1, xsize(2)
            if (istret == 0) y = (j + xstart(2) - 2) * dy
            if (istret /= 0) y = yp(j + xstart(2) - 1)
            um1 = half * (one + tanh((y - inflow_offset) / delta_low))
            um2 = half * (one - tanh((y - inflow_offset - one) / delta_up))
            um = um1 * um2
            ubase = um1 * ((u1 - u2) * um2 + u2)
            bxx1(j, k) = u2 + ramp * (ubase - u2) + (bcx(j, k) - half) * um * inflow_noise * ramp
            bxy1(j, k) = (bcy(j, k) - half) * um * inflow_noise * ramp
            bxz1(j, k) = (bcz(j, k) - half) * um * inflow_noise * ramp
         end do
      end do

      if (iscalar == 1) then
         do is = 1, numscalar
            do k = 1, xsize(3)
               do j = 1, xsize(2)
                  if (istret == 0) y = (j + xstart(2) - 2) * dy
                  if (istret /= 0) y = yp(j + xstart(2) - 1)
                  phi(1, j, k, is) = ramp * cp(is) * half * (tanh((y - inflow_offset) / delta_low) - tanh((y - inflow_offset - one) / delta_up))
               end do
            end do
         end do
      end if

      return
   end subroutine inflow
  !############################################################################

   subroutine outflow(ux, uy, uz, phi)

      use param
      use variables
      use mpi
      use ibm_param

      implicit none

      integer :: i, j, k, code, it, jmax, iloc
      integer :: nouter_local, nouter_global
      integer, dimension(3) :: ix_target
      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux, uy, uz
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi
      real(mytype) :: udx, cx, uxmin, uxmax, uxmean
      real(mytype) :: y, u_inlet_outer_min_local, u_inlet_outer_min
      real(mytype) :: v_top_absmax_local, v_top_absmax
      real(mytype) :: uline_max, ypeak_local, ypeak_global
      real(mytype), dimension(3) :: y_ridge

      udx = one / dx

      uxmax = -1609._mytype
      uxmin = 1609._mytype
      uxmean = zero

      u_inlet_outer_min_local = huge(one)
      v_top_absmax_local = zero
      nouter_local = 0

      ix_target(1) = max(2, int(0.25_mytype * real(nx, mytype)))
      ix_target(2) = max(2, int(0.50_mytype * real(nx, mytype)))
      ix_target(3) = max(2, int(0.75_mytype * real(nx, mytype)))
      y_ridge = -one

      do k = 1, xsize(3)
         do j = 1, xsize(2)
            if (istret == 0) y = (j + xstart(2) - 2) * dy
            if (istret /= 0) y = yp(j + xstart(2) - 1)
            if (ux(nx - 1, j, k) > uxmax) uxmax = ux(nx - 1, j, k)
            if (ux(nx - 1, j, k) < uxmin) uxmin = ux(nx - 1, j, k)
            uxmean = uxmean + ux(nx - 1, j, k)

            if (y > inflow_offset + one + two * delta_up) then
               u_inlet_outer_min_local = min(u_inlet_outer_min_local, ux(1, j, k))
               nouter_local = nouter_local + 1
            end if
         end do
      end do

      do k = 1, xsize(3)
         do i = 1, xsize(1)
            v_top_absmax_local = max(v_top_absmax_local, abs(uy(i, xsize(2) - 1, k)))
         end do
      end do

      do it = 1, 3
         ypeak_local = -one
         if (ix_target(it) >= xstart(1) .and. ix_target(it) <= xend(1)) then
            iloc = ix_target(it) - xstart(1) + 1
            uline_max = -huge(one)
            jmax = 1
            do j = 1, xsize(2)
               if (ux(iloc, j, 1) > uline_max) then
                  uline_max = ux(iloc, j, 1)
                  jmax = j
               end if
            end do
            if (istret == 0) ypeak_local = (jmax + xstart(2) - 2) * dy
            if (istret /= 0) ypeak_local = yp(jmax + xstart(2) - 1)
         end if
         call MPI_ALLREDUCE(ypeak_local, ypeak_global, 1, real_type, MPI_MAX, MPI_COMM_WORLD, code)
         y_ridge(it) = ypeak_global
      end do

      call MPI_ALLREDUCE(MPI_IN_PLACE, uxmax, 1, real_type, MPI_MAX, MPI_COMM_WORLD, code)
      call MPI_ALLREDUCE(MPI_IN_PLACE, uxmin, 1, real_type, MPI_MIN, MPI_COMM_WORLD, code)
      call MPI_ALLREDUCE(MPI_IN_PLACE, uxmean, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
      call MPI_ALLREDUCE(MPI_IN_PLACE, nouter_local, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, code)
      nouter_global = nouter_local
      if (nouter_global > 0) then
         call MPI_ALLREDUCE(MPI_IN_PLACE, u_inlet_outer_min_local, 1, real_type, MPI_MIN, MPI_COMM_WORLD, code)
      else
         u_inlet_outer_min_local = zero
      end if
      call MPI_ALLREDUCE(MPI_IN_PLACE, v_top_absmax_local, 1, real_type, MPI_MAX, MPI_COMM_WORLD, code)

      u_inlet_outer_min = u_inlet_outer_min_local
      v_top_absmax = v_top_absmax_local
      uxmean = uxmean / real(ny * nz, mytype)

      ! Keep the explicit convective BC stable in strongly transient stages.
      cx = max(zero, min(one, uxmax * gdt(itr) * udx))

      do k = 1, xsize(3)
         do j = 1, xsize(2)
            bxxn(j, k) = ux(nx, j, k) - cx * (ux(nx, j, k) - ux(nx - 1, j, k))
            bxyn(j, k) = uy(nx, j, k) - cx * (uy(nx, j, k) - uy(nx - 1, j, k))
            bxzn(j, k) = uz(nx, j, k) - cx * (uz(nx, j, k) - uz(nx - 1, j, k))
         end do
      end do

      if (iscalar == 1) then
         do k = 1, xsize(3)
            do j = 1, xsize(2)
               phi(nx, j, k, :) = phi(nx, j, k, :) - cx * (phi(nx, j, k, :) - phi(nx - 1, j, k, :))
            end do
         end do
      end if

      if (nrank == 0 .and. (mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) &
         write(*, *) "Outflow velocity min max mean=", real(uxmin, 4), real(uxmax, 4), real(uxmean, 4)
      if (nrank == 0 .and. (mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) then
         write(*, *) "DIAG inlet_outer_u_min=", real(u_inlet_outer_min, 4), "top_v_absmax=", real(v_top_absmax, 4)
         write(*, *) "DIAG ridge_y(x=0.25L,0.50L,0.75L)=", real(y_ridge(1), 4), real(y_ridge(2), 4), real(y_ridge(3), 4)
         write(*, *) "DIAG ridge_offset_from_inflow=", real(y_ridge(1)-inflow_offset, 4), real(y_ridge(2)-inflow_offset, 4), real(y_ridge(3)-inflow_offset, 4)
      end if

      return
   end subroutine outflow
  !############################################################################

   subroutine init_walljet(ux1, uy1, uz1, phi1)

      use decomp_2d_io
      use variables
      use param
      use mpi

      implicit none

      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

      integer :: i, j, k, is, code, ierr, ii
      integer, allocatable :: seed(:)
      logical :: has_ux, has_uy, has_uz
      real(mytype) :: x, xnorm, y, ytop_blend, um1, um2, um, ramp0
      real(mytype) :: rux, ruy, ruz
      real(mytype) :: delta_l_local, delta_u_local, u1_local
      real(mytype), parameter :: develop_rate = 0.40_mytype
      real(mytype), parameter :: decay_rate = 0.10_mytype

      ux1 = zero
      uy1 = zero
      uz1 = zero

      if (iread_init == 1) then
         inquire(file = "initial/ux.bin", exist = has_ux)
         inquire(file = "initial/uy.bin", exist = has_uy)
         inquire(file = "initial/uz.bin", exist = has_uz)
         if (.not. (has_ux .and. has_uy .and. has_uz)) then
            if (nrank == 0) then
               write(*, *) 'Error: missing required initial files in ./initial/'
               if (.not. has_ux) write(*, *) '  - ./initial/ux.bin'
               if (.not. has_uy) write(*, *) '  - ./initial/uy.bin'
               if (.not. has_uz) write(*, *) '  - ./initial/uz.bin'
            end if
            call MPI_ABORT(MPI_COMM_WORLD, code, ierr)
            stop
         end if

         if (nrank == 0) write(*, *) 'reading : ', './initial/ux.bin'
         call decomp_2d_read_one(1, ux1, "initial", "ux.bin", io_jet, reduce_prec = .false.)
         if (nrank == 0) write(*, *) 'reading : ', './initial/uy.bin'
         call decomp_2d_read_one(1, uy1, "initial", "uy.bin", io_jet, reduce_prec = .false.)
         if (nrank == 0) write(*, *) 'reading : ', './initial/uz.bin'
         call decomp_2d_read_one(1, uz1, "initial", "uz.bin", io_jet, reduce_prec = .false.)
      else
         ramp0 = one
         if (walljet_ramp_tau > tiny(one)) then
            ramp0 = one - exp(-dt / walljet_ramp_tau)
         end if

         if (iin /= 0 .and. init_noise > zero) then
            call system_clock(count = code)
            if (iin == 2) code = 0
            call random_seed(size = ii)
            allocate(seed(ii))
            do i = 1, ii
               seed(i) = code + 63946 * (nrank + 1) * (i - 1)
            end do
            call random_seed(put = seed)
            deallocate(seed)
            call random_number(ux1)
            call random_number(uy1)
            call random_number(uz1)
         end if

         do k = 1, xsize(3)
            do j = 1, xsize(2)
               if (istret == 0) y = (j + xstart(2) - 2) * dy
               if (istret /= 0) y = yp(j + xstart(2) - 1)
               do i = 1, xsize(1)
                  x = (i + xstart(1) - 2) * dx
                  xnorm = min(one, max(zero, x / max(xlx, tiny(one))))

                  ! Approximate streamwise development to reduce startup transients.
                  delta_l_local = delta_low * (one + develop_rate * xnorm)
                  delta_u_local = delta_up * (one + develop_rate * xnorm)
                  u1_local = max(u2, u1 / (one + decay_rate * xnorm))

                  um1 = half * (one + tanh((y - inflow_offset) / delta_l_local))
                  um2 = half * (one - tanh((y - inflow_offset - one) / delta_u_local))
                  um = um1 * um2
                  ytop_blend = min(one, max(zero, (y - 0.80_mytype * yly) / (0.20_mytype * yly)))

                  rux = ux1(i, j, k)
                  ruy = uy1(i, j, k)
                  ruz = uz1(i, j, k)

                  ux1(i, j, k) = u2 + ramp0 * (um1 * ((u1_local - u2) * um2 + u2) - u2)
                  ux1(i, j, k) = (one - ytop_blend) * ux1(i, j, k) + ytop_blend * u2
                  uy1(i, j, k) = zero
                  uz1(i, j, k) = zero

                  if (iin /= 0 .and. init_noise > zero) then
                     ux1(i, j, k) = ux1(i, j, k) + init_noise * um * (two * rux - one)
                     uy1(i, j, k) = init_noise * um * (two * ruy - one)
                     uz1(i, j, k) = init_noise * um * (two * ruz - one)
                  end if
               end do
            end do
         end do

         if (nz == 1) uz1 = zero
      end if

      if (iscalar == 1) then
         do is = 1, numscalar
            phi1(:, :, :, is) = zero
         end do
      end if

      if (nrank == 0) write(*, *) '# init end ok'

      return
   end subroutine init_walljet
  !############################################################################
  
  subroutine momentum_forcing_walljet(dux1, duy1, duz1, ux1, uy1, uz1)

      use param

      implicit none

      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1

      integer :: i, j, k
      real(mytype) :: x, y, ramp, lambda_top, lambda_backflow
      real(mytype), parameter :: y_top_start_ratio = 0.75_mytype
      real(mytype), parameter :: top_damp = 1.20_mytype
      real(mytype), parameter :: backflow_damp = 2.50_mytype

      ramp = one
      if (walljet_ramp_tau > tiny(one)) then
         ramp = one - exp(-max(zero, real(itime - 1, mytype) * dt) / walljet_ramp_tau)
      end if

      do k = 1, xsize(3)
         do j = 1, xsize(2)
            if (istret == 0) y = (j + xstart(2) - 2) * dy
            if (istret /= 0) y = yp(j + xstart(2) - 1)

            ! Damping layer near the top boundary to absorb large coherent startup vortices.
            lambda_top = zero
            if (y > y_top_start_ratio * yly) then
               lambda_top = (y - y_top_start_ratio * yly) / ((one - y_top_start_ratio) * yly)
               lambda_top = max(zero, min(one, lambda_top))
               lambda_top = top_damp * lambda_top * lambda_top
            end if

            do i = 1, xsize(1)
               x = (i + xstart(1) - 2) * dx

               if (lambda_top > zero) then
                  dux1(i, j, k, 1) = dux1(i, j, k, 1) + ramp * lambda_top * (u2 - ux1(i, j, k))
                  duy1(i, j, k, 1) = duy1(i, j, k, 1) + ramp * lambda_top * (zero - uy1(i, j, k))
                  duz1(i, j, k, 1) = duz1(i, j, k, 1) + ramp * lambda_top * (zero - uz1(i, j, k))
               end if

               ! Suppress large-scale recirculating backflow above the jet near inlet.
               if (x < 0.35_mytype * xlx .and. y > inflow_offset + one + two * delta_up) then
                  if (ux1(i, j, k) < zero) then
                     lambda_backflow = backflow_damp * (one - x / max(0.35_mytype * xlx, tiny(one)))
                     lambda_backflow = max(zero, lambda_backflow)
                     dux1(i, j, k, 1) = dux1(i, j, k, 1) - ramp * lambda_backflow * ux1(i, j, k)
                  end if
               end if
            end do
         end do
      end do

  end subroutine momentum_forcing_walljet

  !############################################################################
  subroutine postprocess_walljet(ux1, uy1, uz1, pp3, phi1, ep1)

    use var, ONLY : nzmsize
    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

  end subroutine postprocess_walljet

  !############################################################################
  subroutine visu_walljet_init (visu_initialised)

    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "critq", 1, 0, output2D, mytype)

    visu_initialised = .true.
    
  end subroutine visu_walljet_init
  
  !############################################################################
  subroutine visu_walljet(ux1, uy1, uz1, pp3, phi1, ep1, num)

    use var, only : ux2, uy2, uz2, ux3, uy3, uz3
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    use var, ONLY : nxmsize, nymsize, nzmsize
    use visu, only : write_field
    use ibm_param, only : ubcx,ubcy,ubcz

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    integer, intent(in) :: num

  end subroutine visu_walljet
  !############################################################################
  
end module walljet
