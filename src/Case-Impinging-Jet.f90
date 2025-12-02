!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module impingingjet

  USE decomp_2d_constants
  USE decomp_2d_mpi
  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  character(len=*), parameter :: io_jet = "io-jet"

  PRIVATE !! All functions/subroutines private by default
  PUBLIC :: boundary_conditions_impingingjet, momentum_forcing_impingingjet, init_impingingjet, &
            postprocess_impingingjet, visu_impingingjet, visu_impingingjet_init

  ! T. Dairay, V. Fortuné, E. Lamballais, and L.-E. Brizzi, 
  ! Direct numerical simulation of a turbulent jet impinging on a heated wall, 
  ! J. Fluid Mech. 764, 362 (2015).

contains

  subroutine boundary_conditions_impingingjet (ux, uy, uz, phi)
    USE param
    USE variables

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux, uy, uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    integer :: is, isT
    real(mytype) :: qflux

    !Top boundary (yn)
    call inflow (ux, uy, uz)

    !Side boundary (x1,xn,z1,zn)
    call outflow (ux, uy, uz, phi)

    !Scalar boundary conditions
    !Top boundary (yn)
    if (iscalar.eq.1) then
       if (xend(2).eq.ny) then
          !set all scalar components at top to 0
          do is = 1, numscalar
             phi(:,xsize(2),:,is) = zero
          enddo
       endif
    endif

    !Bottom boundary (y1)
    if (iscalar.eq.1) then
       if (xstart(2).eq.1) then
          isT = 1; qflux = one-exp(-half*t*t)
          call scalar_neumann (phi,isT,qflux)
          do is = 2, numscalar
             phi(:,xstart(2),:,is) = one-exp(-half*t*t)
          enddo
       endif
    endif

    return
  end subroutine boundary_conditions_impingingjet
  !########################################################################

  subroutine inflow (ux, uy, uz)
    USE param
    USE variables

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux, uy, uz

    integer :: i, k, is, m
    real(mytype) :: x, z, r, theta, delta, sw
    real(mytype) :: um, sigma, Stc, omega_c, tauA, sigmaA, rhoA, xi
    real(mytype) :: uxprime, uyprime, uzprime

    integer, parameter :: NMODE = 8
    real(mytype), save :: Amx(NMODE), Amy(NMODE), Amz(NMODE)
    real(mytype), save :: phx(NMODE), phy(NMODE), phz(NMODE)
    real(mytype), save :: omegax(NMODE), omegay(NMODE), omegaz(NMODE)
    logical, save      :: inflow_init = .false.
    
    delta = three*sqrt(dx*dz)

    sigma = 100._mytype
    Stc = 1._mytype
    omega_c = two*pi*Stc
    tauA = one/omega_c
    rhoA = exp(-dt/tauA)

    sigmaA = inflow_noise*sqrt(two/real(NMODE,mytype))
    
    if (.not. inflow_init) then
      call random_seed()

      call random_number(Amx);  Amx = two*Amx - one  ! map to [-1,1]
      call random_number(Amy);  Amy = two*Amy - one
      call random_number(Amz);  Amz = two*Amz - one

      call random_number(phx);  phx = phx * two*pi
      call random_number(phy);  phy = phy * two*pi
      call random_number(phz);  phz = phz * two*pi

      call random_number(omegax); omegax = (two*omegax - one) * omega_c
      call random_number(omegay); omegay = (two*omegay - one) * omega_c
      call random_number(omegaz); omegaz = (two*omegaz - one) * omega_c
 
      inflow_init = .true.
    endif

    if (iin.ne.0) then
      ! if (nrank .eq. 0) write(*,*) "# inflow using synthetic perturbation: A*cos(m*theta+phi)"
      do m = 1, NMODE
         call randn_gauss(xi)
         Amx(m) = rhoA*Amx(m)+sqrt(one-rhoA*rhoA)*sigmaA*xi
         call randn_gauss(xi)
         Amy(m) = rhoA*Amy(m)+sqrt(one-rhoA*rhoA)*sigmaA*xi
         call randn_gauss(xi)
         Amz(m) = rhoA*Amz(m)+sqrt(one-rhoA*rhoA)*sigmaA*xi
 
         phx(m) = modulo(phx(m)+omegax(m)*dt,two*pi)
         phy(m) = modulo(phy(m)+omegay(m)*dt,two*pi)
         phz(m) = modulo(phz(m)+omegaz(m)*dt,two*pi)
      end do
    endif

    !Top boundary
    do k = 1,xsize(3)
       z = (k+xstart(3)-2)*dz-zlz/two
       do i = 1,xsize(1)
          x = (i+xstart(1)-2)*dx-xlx/two
          r = sqrt(x*x+z*z)
          theta = atan2(z,x)

          uxprime = zero; uyprime = zero; uzprime = zero
          if (iin.ne.0) then
            do m = 1, NMODE
               uxprime = uxprime + Amx(m)*cos(m*theta+phx(m))
               uyprime = uyprime + Amy(m)*cos(m*theta+phy(m))
               uzprime = uzprime + Amz(m)*cos(m*theta+phz(m))
            enddo
          endif

          if (r.le.half) then
             um = 0.7_mytype*sigma*exp(one)*(half-r)*exp(-sigma*(half-r))
             sw = half*(one+cos(pi*(r-(half-delta))/delta))
             if (r.le.half-delta) sw = one
             byxn(i,k) = um*uxprime*sw
             byyn(i,k) = -(u2+two)/u2*(one-(two*r)**u2)*sw + um*uyprime*sw
             byzn(i,k) = um*uzprime*sw
          else
             byxn(i,k) = zero
             byyn(i,k) = zero
             byzn(i,k) = zero
          endif
       enddo
    enddo

    return

  contains

    subroutine randn_gauss(xx)
       real(mytype), intent(out) :: xx
       real(mytype) :: r1, r2
       call random_number(r1)
       call random_number(r2)
       xx = sqrt(-two*log(r1))*cos(two*pi*r2)
    end subroutine randn_gauss

  end subroutine inflow
  !########################################################################

  subroutine outflow (ux, uy, uz, phi)
    USE param
    USE variables
    USE MPI
    USE ibm_param

    implicit none

    integer :: i, j, k, is, ierr
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux, uy, uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype) :: udx, udy, udz, uddx, uddy, uddz
    real(mytype) :: x, y, z, r, un, cx1, cxn, cz1, czn
    real(mytype) :: ux1min, ux1max, ux1mean, uxnmin, uxnmax, uxnmean
    real(mytype) :: uz1min, uz1max, uz1mean, uznmin, uznmax, uznmean

    udx=one/dx; udy=one/dy; udz=one/dz; uddx=half/dx; uddy=half/dy; uddz=half/dz

    ux1max = -1609._mytype; ux1min = 1609._mytype; ux1mean = zero
    uxnmax = -1609._mytype; uxnmin = 1609._mytype; uxnmean = zero
    uz1max = -1609._mytype; uz1min = 1609._mytype; uz1mean = zero
    uznmax = -1609._mytype; uznmin = 1609._mytype; uznmean = zero

    !x1 (-ex)
    do k = 1, xsize(3)
       do j = 1, xsize(2)
          if (ux(2,j,k).gt.ux1max) ux1max = ux(2,j,k)
          if (ux(2,j,k).lt.ux1min) ux1min = ux(2,j,k)
          ux1mean = ux1mean + ux(2,j,k)
       enddo
    enddo
    ux1mean = ux1mean/xsize(2)/xsize(3)
    call MPI_ALLREDUCE(MPI_IN_PLACE,ux1max,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,ux1min,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,ux1mean,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    ux1mean = ux1mean/nproc
    if (iopen.eq.0) then
       !convective outflow using parabolic profile
       do k = 1, xsize(3)
          z = (k+xstart(3)-2)*dz-zlz/two
          do j = 1, xsize(2)
             if (istret.eq.0) y=(j+xstart(2)-2)*dy-yly/two
             if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/two
             x = -xlx/two
             r = sqrt(x*x+z*z)
             un = three*(one-four*y*y/yly/yly)/yly/r/sixteen*x/r
             cx1 = un*gdt(itr)*udx
             bxx1(j,k) = ux(1,j,k)-cx1*(ux(2,j,k)-ux(1,j,k))
             bxy1(j,k) = uy(1,j,k)-cx1*(uy(2,j,k)-uy(1,j,k))
             bxz1(j,k) = uz(1,j,k)-cx1*(uz(2,j,k)-uz(1,j,k))
             if (iscalar.eq.1) then
                phi(1,j,k,:) = phi(1,j,k,:)-cx1*(phi(2,j,k,:)-phi(1,j,k,:))
             endif
          enddo
       enddo
    elseif (iopen.eq.1) then
       !open boundary: zero-gradient + reverse flow (inflow velocity set to 0)
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             un = -ux(2,j,k)
             if (un.ge.zero) then
                bxx1(j,k) = ux(2,j,k)
                bxy1(j,k) = uy(2,j,k)
                bxz1(j,k) = uz(2,j,k)
                if (iscalar.eq.1) then
                   phi(1,j,k,:) = phi(2,j,k,:)
                endif
             else
                bxx1(j,k) = zero
                bxy1(j,k) = zero
                bxz1(j,k) = zero
                if (iscalar.eq.1) then
                   phi(1,j,k,:) = zero
                endif
             endif
          enddo
       enddo
    endif

    !xn (+ex)
    do k = 1, xsize(3)
       do j = 1, xsize(2)
          if (ux(nx-1,j,k).gt.uxnmax) uxnmax = ux(nx-1,j,k)
          if (ux(nx-1,j,k).lt.uxnmin) uxnmin = ux(nx-1,j,k)
          uxnmean = uxnmean + ux(nx-1,j,k)
       enddo
    enddo
    uxnmean = uxnmean/xsize(2)/xsize(3)
    call MPI_ALLREDUCE(MPI_IN_PLACE,uxnmax,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,uxnmin,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,uxnmean,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    uxnmean = uxnmean/nproc
    if (iopen.eq.0) then
       !convective outflow using parabolic profile
       do k = 1, xsize(3)
          z = (k+xstart(3)-2)*dz-zlz/two
          do j = 1, xsize(2)
             if (istret.eq.0) y=(j+xstart(2)-2)*dy-yly/two
             if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/two
             x = xlx/two
             r = sqrt(x*x+z*z)
             un = three*(one-four*y*y/yly/yly)/yly/r/sixteen*x/r
             cxn = un*gdt(itr)*udx
             bxxn(j,k) = ux(nx,j,k)-cxn*(ux(nx,j,k)-ux(nx-1,j,k))
             bxyn(j,k) = uy(nx,j,k)-cxn*(uy(nx,j,k)-uy(nx-1,j,k))
             bxzn(j,k) = uz(nx,j,k)-cxn*(uz(nx,j,k)-uz(nx-1,j,k))
             if (iscalar.eq.1) then
                phi(nx,j,k,:) = phi(nx,j,k,:)-cxn*(phi(nx,j,k,:)-phi(nx-1,j,k,:))
             endif
          enddo
       enddo
    elseif (iopen.eq.1) then
       !open boundary: zero-gradient + reverse flow (inflow velocity set to 0)
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             un = ux(nx-1,j,k)
             if (un.ge.zero) then
                bxxn(j,k) = ux(nx-1,j,k)
                bxyn(j,k) = uy(nx-1,j,k)
                bxzn(j,k) = uz(nx-1,j,k)
                if (iscalar.eq.1) then
                   phi(nx,j,k,:) = phi(nx-1,j,k,:)
                endif
             else
                bxxn(j,k) = zero
                bxyn(j,k) = zero
                bxzn(j,k) = zero
                if (iscalar.eq.1) then
                   phi(nx,j,k,:) = zero
                endif
             endif
          enddo
       enddo
    endif

    !z1 (-ez)
    if (xstart(3).eq.1) then
       do j = 1, xsize(2)
          do i = 1, xsize(1)
             if (uz(i,j,2).gt.uz1max) uz1max = uz(i,j,2)
             if (uz(i,j,2).lt.uz1min) uz1min = uz(i,j,2)
             uz1mean = uz1mean + uz(i,j,2)
          enddo
       enddo
    endif
    uz1mean = uz1mean/xsize(1)/xsize(2)
    call MPI_ALLREDUCE(MPI_IN_PLACE,uz1max,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,uz1min,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,uz1mean,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    uz1mean = uz1mean/nproc
    if (xstart(3).eq.1) then
       if (iopen.eq.0) then
          !convective outflow using parabolic profile
          do j = 1, xsize(2)
             if (istret.eq.0) y=(j+xstart(2)-2)*dy-yly/two
             if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/two
             do i = 1, xsize(1)
                x = (i+xstart(1)-2)*dx-xlx/two
                z = -zlz/two
                r = sqrt(x*x+z*z)
                un = three*(one-four*y*y/yly/yly)/yly/r/sixteen*z/r
                cz1 = un*gdt(itr)*udz
                bzx1(i,j) = ux(i,j,1)-cz1*(ux(i,j,2)-ux(i,j,1))
                bzy1(i,j) = uy(i,j,1)-cz1*(uy(i,j,2)-uy(i,j,1))
                bzz1(i,j) = uz(i,j,1)-cz1*(uz(i,j,2)-uz(i,j,1))
                if (iscalar.eq.1) then
                   phi(i,j,1,:) = phi(i,j,1,:)-cz1*(phi(i,j,2,:)-phi(i,j,1,:))
                endif
             enddo
          enddo
       elseif (iopen.eq.1) then
          !open boundary: zero-gradient + reverse flow (inflow velocity set to 0)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                un = -uz(i,j,2)
                if (un.ge.zero) then
                   bzx1(i,j) = ux(i,j,2)
                   bzy1(i,j) = uy(i,j,2)
                   bzz1(i,j) = uz(i,j,2)
                   if (iscalar.eq.1) then
                      phi(i,j,1,:) = phi(i,j,2,:)
                   endif
                else
                   bzx1(i,j) = zero
                   bzy1(i,j) = zero
                   bzz1(i,j) = zero
                   if (iscalar.eq.1) then
                      phi(i,j,1,:) = zero
                   endif
                endif
             enddo
          enddo
       endif
    endif

    !zn (+ez)
    if (xend(3).eq.nz) then
       do j = 1, xsize(2)
          do i = 1, xsize(1)
             if (uz(i,j,xsize(3)-1).gt.uznmax) uznmax = uz(i,j,xsize(3)-1)
             if (uz(i,j,xsize(3)-1).lt.uznmin) uznmin = uz(i,j,xsize(3)-1)
             uznmean = uznmean + uz(i,j,xsize(3)-1)
          enddo
       enddo
    endif
    uznmean = uznmean/xsize(1)/xsize(2)
    call MPI_ALLREDUCE(MPI_IN_PLACE,uznmax,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,uznmin,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,uznmean,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    uznmean = uznmean/nproc
    if (xend(3).eq.nz) then
       if (iopen.eq.0) then
          !convective outflow using parabolic profile
          do j = 1, xsize(2)
             if (istret.eq.0) y=(j+xstart(2)-2)*dy-yly/two
             if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/two
             do i = 1, xsize(1)
                x = (i+xstart(1)-2)*dx-xlx/two
                z = zlz/two
                r = sqrt(x*x+z*z)
                un = three*(one-four*y*y/yly/yly)/yly/r/sixteen*z/r
                czn = un*gdt(itr)*udz
                bzxn(i,j) = ux(i,j,nz)-czn*(ux(i,j,nz)-ux(i,j,nz-1))
                bzyn(i,j) = uy(i,j,nz)-czn*(uy(i,j,nz)-uy(i,j,nz-1))
                bzzn(i,j) = uz(i,j,nz)-czn*(uz(i,j,nz)-uz(i,j,nz-1))
                if (iscalar.eq.1) then
                   phi(i,j,nz,:) = phi(i,j,nz,:)-czn*(phi(i,j,nz,:)-phi(i,j,nz-1,:))
                endif
             enddo
          enddo
       elseif (iopen.eq.1) then
          !open boundary: zero-gradient + reverse flow (inflow velocity set to 0)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                un = uz(i,j,xsize(3)-1)
                if (un.ge.zero) then
                   bzxn(i,j) = ux(i,j,xsize(3)-1)
                   bzyn(i,j) = uy(i,j,xsize(3)-1)
                   bzzn(i,j) = uz(i,j,xsize(3)-1)
                   if (iscalar.eq.1) then
                      phi(i,j,nz,:) = phi(i,j,nz-1,:)
                   endif
                else
                   bzxn(i,j) = zero
                   bzyn(i,j) = zero
                   bzzn(i,j) = zero
                   if (iscalar.eq.1) then
                      phi(i,j,nz,:) = zero
                   endif
                endif
             enddo
          enddo
       endif
    endif

    if (nrank==0 .and. (mod(itime, ilist)==0 .or. itime==ifirst .or. itime==ilast)) then
       if (iopen.eq.0) then
          write(*,*) "Convective outflow BCs: using parabolic profile"
       elseif (iopen.eq.1) then
          write(*,*) "Open boundary BCs: zero-gradient + reverse flow"
       endif
       write(*,*) "Outflow velocity ux1 min max=",real(ux1min,4),real(ux1max,4),real(ux1mean,4)
       write(*,*) "Outflow velocity uxn min max=",real(uxnmin,4),real(uxnmax,4),real(uxnmean,4)
       write(*,*) "Outflow velocity uz1 min max=",real(uz1min,4),real(uz1max,4),real(uz1mean,4)
       write(*,*) "Outflow velocity uzn min max=",real(uznmin,4),real(uznmax,4),real(uznmean,4)
    endif

    return
  end subroutine outflow
  !########################################################################

  subroutine scalar_neumann (phi,is,qflux)
    USE param
    USE variables
 
    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    integer,intent(in) :: is
    real(mytype),intent(in):: qflux

    integer :: i, k
    real(mytype) :: dy0, dy1, a0, a1, a2

    !Dairay et al., Int. J. Heat Fluid Flow (2014) 
    if (xstart(2).eq.1) then
       if (istret.eq.0) then
          dy0 = dy
          do k = 1, xsize(3)
             do i = 1, xsize(1)
                !T1 = (18*T2 - 9*T3 + 2*T4 + 6*dy0*qflux)/11
                phi(i,1,k,is) = (eighteen*phi(i,2,k,is)-nine*phi(i,3,k,is)+two*phi(i,4,k,is)+six*dy0*qflux)/eleven
             enddo
          enddo
          if (nrank .eq. 0) write(*,*) 'Neumann BC for scalar (third-order uniform grid)'
       else
          dy0 = yp(xstart(2)+1)-yp(xstart(2))
          dy1 = yp(xstart(2)+2)-yp(xstart(2)+1)
          a0 = -(two*dy0 + dy1)/dy0/(dy0 + dy1)
          a1 =  (dy0+dy1)/dy0/dy1
          a2 = - dy0/dy1/(dy0 + dy1)
          do k = 1, xsize(3)
             do i = 1, xsize(1)
                !T1 = (-qflux - a1*T2 - a2*T3)/a0
                phi(i,1,k,is) = -(qflux+a1*phi(i,2,k,is)+a2*phi(i,3,k,is))/a0
             enddo
          enddo
          if (nrank .eq. 0) write(*,*) 'Neumann BC for scalar (second-order stretched grid)'
       endif
    endif

  end subroutine scalar_neumann
  !########################################################################

  subroutine momentum_forcing_impingingjet(dux1,duy1,duz1,ux1,uy1,uz1)
    USE param
    USE variables

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3), ntime) :: dux1, duy1, duz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1

    integer :: i, j, k
    real(mytype) :: x, y, z, r, um, lambda

    if (ifringe.ne.0) then
       do k = 1,xsize(3)
          z = (k+xstart(3)-2)*dz-zlz/two
          do j = 1,xsize(2)
             if (istret.eq.0) y=(j+xstart(2)-2)*dy-yly/two
             if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/two
             do i = 1,xsize(1)
                x = (i+xstart(1)-2)*dx-xlx/two
                r = sqrt(x*x+z*z)
                if (r.ge.fringe_rm) then
                   lambda = half*(one+tanh(fringe_beta*(r-fringe_rm)-four))
                   um = three*(one-four*y*y/yly/yly)/yly/r/sixteen
                   if (ifringe.ne.2) then
                      duy1(i,j,k,1) = duy1(i,j,k,1) + lambda*(zero-uy1(i,j,k))
                   endif
                   !ifringe=2 means only damp u and w
                   dux1(i,j,k,1) = dux1(i,j,k,1) + lambda*(um*x/r-ux1(i,j,k))
                   duz1(i,j,k,1) = duz1(i,j,k,1) + lambda*(um*z/r-uz1(i,j,k))
                endif
             enddo
          enddo
       enddo
      !  if (nrank .eq. 0) write(*,*) '# fringe forcing applied at rm=', fringe_rm
    endif

    return
  end subroutine momentum_forcing_impingingjet
  !########################################################################

  subroutine init_impingingjet (ux1,uy1,uz1,phi1)
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    integer :: is

    ux1=zero; uy1=zero; uz1=zero

    if(iin.eq.4) then
      !Read velocity field
      if (nrank.eq.0) write(*,*) 'reading : ', './data/ux.bin'
      call decomp_2d_read_one(1,ux1,'data','ux.bin',io_jet,reduce_prec=.false.)
      if (nrank.eq.0) write(*,*) 'reading : ', './data/uy.bin'
      call decomp_2d_read_one(1,uy1,'data','uy.bin',io_jet,reduce_prec=.false.)
      if (nrank.eq.0) write(*,*) 'reading : ', './data/uz.bin'
      call decomp_2d_read_one(1,uz1,'data','uz.bin',io_jet,reduce_prec=.false.)
    endif

    if (iscalar.eq.1) then
      do is = 1, numscalar
          phi1(:,:,:,is) = zero
      enddo
    endif
    
    if (nrank .eq. 0) write(*,*) '# init end ok'

    return
  end subroutine init_impingingjet
  !########################################################################

  subroutine postprocess_impingingjet (ux1,uy1,uz1,pp3,phi1,ep1)
    use var, ONLY : nzmsize
    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

  end subroutine postprocess_impingingjet

  subroutine visu_impingingjet_init (visu_initialised)

    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "critq", 1, 0, output2D, mytype)

    visu_initialised = .true.
    
  end subroutine visu_impingingjet_init
  !########################################################################

  subroutine visu_impingingjet (ux1, uy1, uz1, pp3, phi1, ep1, num)
    use var, only : ux2, uy2, uz2, ux3, uy3, uz3
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2
    USE var, only : ta3,tb3,tc3,td3,te3,tf3,di3
    use var, ONLY : nxmsize, nymsize, nzmsize
    use visu, only : write_field, write_xdmf_vector
    use ibm_param, only : ubcx,ubcy,ubcz

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    integer, intent(in) :: num

    integer :: i, k, ierr
    real(mytype) :: dy0, dy1, numax01, numax02
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: nu01, nu02

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
    di1(:,:,:) = tf1(:,:,:)-th1(:,:,:)
    call write_field(di1, ".", "vorx", num, flush = .true.)
    di1 = zero
    di1(:,:,:) = tg1(:,:,:)-tc1(:,:,:)
    call write_field(di1, ".", "vory", num, flush = .true.)
    di1 = zero
    di1(:,:,:) = tb1(:,:,:)-td1(:,:,:)
    call write_field(di1, ".", "vorz", num, flush = .true.)
    call write_xdmf_vector(".", "Vort", "vorx", "vory", "vorz", num)
   !  di1(:,:,:)=sqrt(  (tf1(:,:,:)-th1(:,:,:))**2 &
   !                  + (tg1(:,:,:)-tc1(:,:,:))**2 &
   !                  + (tb1(:,:,:)-td1(:,:,:))**2)
   !  call write_field(di1, ".", "vort", num, flush = .true.) ! Reusing temporary array, force flush

    !Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
    di1 = zero
    di1(:,:,:) = - half*(ta1(:,:,:)**2 + te1(:,:,:)**2 + ti1(:,:,:)**2) &
                 - td1(:,:,:) * tb1(:,:,:) &
                 - tg1(:,:,:) * tc1(:,:,:) &
                 - th1(:,:,:) * tf1(:,:,:)
    call write_field(di1, ".", "critq", num, flush = .true.) ! Reusing temporary array, force flush

    !Nusselt number
    numax01 = -1609._mytype; numax02 = -1609._mytype;
    if (iscalar.eq.1) then
       nu01 = zero; nu02 = zero
       if (xstart(2).eq.1) then
          do k = 1, xsize(3)
             do i = 1, xsize(1)
                if (abs(phi1(i,1,k,1)).gt.1.0e-12_mytype) nu01(i,1,k) = one/phi1(i,1,k,1)
                if (nu01(i,1,k).gt.numax01) numax01 = nu01(i,1,k)

                if (numscalar.ge.2) then
                   if (istret.eq.0) then
                      dy0 = dy
                      nu02(i,1,k) = -(-eleven*phi1(i,1,k,2)+eighteen*phi1(i,2,k,2)-nine*phi1(i,3,k,2)+two*phi1(i,4,k,2))/six/dy0
                   else
                      dy0 = yp(xstart(2)+1)-yp(xstart(2))
                      dy1 = yp(xstart(2)+2)-yp(xstart(2)+1)
                      nu02(i,1,k) = -(-phi1(i,1,k,2)*(two*dy0+dy1)/dy0/(dy0+dy1) &
                                    + phi1(i,2,k,2)*(dy0+dy1)/dy0/dy1 - phi1(i,3,k,2)*dy0/dy1/(dy0+dy1))
                   endif
                   if (nu02(i,1,k).gt.numax02) numax02 = nu02(i,1,k)
                endif
             enddo
          enddo
       endif
       
       call write_field(nu01, ".", "nu01", num)
       call MPI_ALLREDUCE(MPI_IN_PLACE,numax01,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
       if (nrank==0) write(*,*) "Nusselt number nu01 max=", real(numax01,4)

       if (numscalar.ge.2) then
          call write_field(nu02, ".", "nu02", num)
          call MPI_ALLREDUCE(MPI_IN_PLACE,numax02,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
          if (nrank==0) write(*,*) "Nusselt number nu02 max=", real(numax02,4)
       endif
    endif

  end subroutine visu_impingingjet
  !########################################################################

end module impingingjet
