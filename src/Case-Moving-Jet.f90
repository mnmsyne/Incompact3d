!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module movingjet

  USE decomp_2d_constants
  USE decomp_2d_mpi
  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  character(len=*), parameter :: io_jet = "io-jet"

  PRIVATE !! All functions/subroutines private by default
  PUBLIC :: boundary_conditions_movingjet, momentum_forcing_movingjet, init_movingjet, &
            postprocess_movingjet, visu_movingjet, visu_movingjet_init

contains

  subroutine boundary_conditions_movingjet (ux, uy, uz, phi)
    USE param
    USE variables

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux, uy, uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    integer :: i, j, k, is

    !Bottom boundary (y1)
    if (iscalar.eq.1) then
       if (xstart(2).eq.1) then
          do is = 1, numscalar
             do k = 1, xsize(3)
                do i = 1, xsize(1)
                   phi(i,xstart(2),k,is) = phi(i,xstart(2)+1,k,is)
                enddo
             enddo
          enddo
       endif
    endif

    !Top boundary (yn)
    call inflow (ux, uy, uz, phi)

    !Side boundary (x1,xn)
    call outflow (ux, uy, uz, phi)

    return
  end subroutine boundary_conditions_movingjet
  !########################################################################

  subroutine inflow (ux, uy, uz, phi)
    USE param
    USE variables

    implicit none

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux, uy, uz
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    integer :: i, j, k, is
    real(mytype) :: x, z, r, um
    real(mytype), dimension(xsize(1),xsize(3)) :: uxprime, uyprime, uzprime

    call random_number(uxprime)
    call random_number(uyprime)
    call random_number(uzprime)

    !Top boundary
    do k = 1,xsize(3)
       z = (k+xstart(3)-2)*dz-zlz/two
       do i = 1,xsize(1)
          x = (i+xstart(1)-2)*dx-xlx/two
          r = sqrt(x*x)
          um = one-(two*r)**u2
          if (r.le.half) then
             byxn(i,k) = zero + (uxprime(i,k)-half)*um*inflow_noise
             byyn(i,k) = -(u2+one)/u2*(one-(two*r)**u2)*(one-exp(-half*t*t)) + (uyprime(i,k)-half)*um*inflow_noise
             byzn(i,k) = zero + (uzprime(i,k)-half)*um*inflow_noise
          else
             byxn(i,k) = zero
             byyn(i,k) = zero
             byzn(i,k) = zero
          endif
       enddo
    enddo

    if (iscalar.eq.1) then
       if (xend(2).eq.ny) then
          do is = 1, numscalar
             do k = 1,xsize(3)
                z = (k+xstart(3)-2)*dz-zlz/two
                do i = 1,xsize(1)
                   x = (i+xstart(1)-2)*dx-xlx/two
                   r = sqrt(x*x)
                   if (r.le.half) then
                      phi(i,xsize(2),k,is) = cp(is)*(one-exp(-half*t*t))
                   else
                      phi(i,xsize(2),k,is) = phi(i,xsize(2)-1,k,is)
                   endif
                enddo
             enddo
          enddo
       endif
    endif

    return
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
    real(mytype) :: y, un, cx1, cxn
    real(mytype) :: ux1min, ux1max, ux1mean, uxnmin, uxnmax, uxnmean

    udx=one/dx; udy=one/dy; udz=one/dz; uddx=half/dx; uddy=half/dy; uddz=half/dz

    ux1max = -1609._mytype; ux1min = 1609._mytype; ux1mean = zero
    uxnmax = -1609._mytype; uxnmin = 1609._mytype; uxnmean = zero

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
          do j = 1, xsize(2)
             if (istret.eq.0) y=(j+xstart(2)-2)*dy-yly/two
             if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/two
             un = -three/four*(one-four*y*y/yly/yly)/yly
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
          do j = 1, xsize(2)
             if (istret.eq.0) y=(j+xstart(2)-2)*dy-yly/two
             if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/two
             un = three/four*(one-four*y*y/yly/yly)/yly
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

    if (nrank==0 .and. (mod(itime, ilist)==0 .or. itime==ifirst .or. itime==ilast)) then
       if (iopen.eq.0) then
          write(*,*) "Convective outflow BCs: using parabolic profile"
       elseif (iopen.eq.1) then
          write(*,*) "Open boundary BCs: zero-gradient + reverse flow"
       endif
       write(*,*) "Outflow velocity ux1 min max mean=",real(ux1min,4),real(ux1max,4),real(ux1mean,4)
       write(*,*) "Outflow velocity uxn min max mean=",real(uxnmin,4),real(uxnmax,4),real(uxnmean,4)
    endif

    return
  end subroutine outflow
  !########################################################################

  subroutine momentum_forcing_movingjet(dux1,duy1,duz1,ux1,uy1,uz1)
    USE param
    USE variables

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3), ntime) :: dux1, duy1, duz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1

    integer :: i, j, k
    real(mytype) :: x, y, um, lambda

    if (ifringe.ne.0) then
       do k = 1,xsize(3)
          do j = 1,xsize(2)
             if (istret.eq.0) y=(j+xstart(2)-2)*dy-yly/two
             if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/two
             do i = 1,xsize(1)
                x = (i+xstart(1)-2)*dx-xlx/two
                if (x.ge.fringe_rm) then
                   lambda = half*(one+tanh(fringe_beta*(x-fringe_rm)-four))
                   um = three/four*(one-four*y*y/yly/yly)/yly
                   dux1(i,j,k,1) = dux1(i,j,k,1) + lambda*(um-ux1(i,j,k))
                   duy1(i,j,k,1) = duy1(i,j,k,1) + lambda*(zero-uy1(i,j,k))
                   duz1(i,j,k,1) = duz1(i,j,k,1) + lambda*(zero-uz1(i,j,k))
                elseif (x.le.-fringe_rm) then
                   lambda = half*(one+tanh(-fringe_beta*(x+fringe_rm)-four))
                   um = -three/four*(one-four*y*y/yly/yly)/yly
                   dux1(i,j,k,1) = dux1(i,j,k,1) + lambda*(um-ux1(i,j,k))
                   duy1(i,j,k,1) = duy1(i,j,k,1) + lambda*(zero-uy1(i,j,k))
                   duz1(i,j,k,1) = duz1(i,j,k,1) + lambda*(zero-uz1(i,j,k))
                endif
             enddo
          enddo
       enddo
    endif

    return
  end subroutine momentum_forcing_movingjet
  !########################################################################

  subroutine init_movingjet (ux1,uy1,uz1,phi1)
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
  end subroutine init_movingjet
  !########################################################################

  subroutine postprocess_movingjet (ux1,uy1,uz1,pp3,phi1,ep1)
    use var, ONLY : nzmsize
    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

  end subroutine postprocess_movingjet

  subroutine visu_movingjet_init (visu_initialised)

    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "critq", 1, 0, output2D, mytype)

    visu_initialised = .true.
    
  end subroutine visu_movingjet_init
  !########################################################################

  subroutine visu_movingjet (ux1, uy1, uz1, pp3, phi1, ep1, num)
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

    !Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
    di1 = zero
    di1(:,:,:) = - half*(ta1(:,:,:)**2 + te1(:,:,:)**2 + ti1(:,:,:)**2) &
                 - td1(:,:,:) * tb1(:,:,:) &
                 - tg1(:,:,:) * tc1(:,:,:) &
                 - th1(:,:,:) * tf1(:,:,:)
    call write_field(di1, ".", "critq", num, flush = .true.) ! Reusing temporary array, force flush

  end subroutine visu_movingjet
  !########################################################################

end module movingjet
