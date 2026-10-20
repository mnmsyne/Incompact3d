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
            visu_swirljet, visu_swirljet_init

contains

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

    integer  :: j,k,is
    real(mytype) :: y, z, r, theta
    real(mytype) :: ratio, um, rm, delta
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    call validate_swirljet_radii()

    delta = 0.1_mytype
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
          if (r.gt.rao) then
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
    endif

    return
  end subroutine inflow

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

    if (u1 == zero) then
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
       if (u2==zero) then
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

    call validate_swirljet_radii()

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
             if (r.gt.rao) then
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

#ifdef DEBG
    if (nrank .eq. 0) write(*,*) '# init end ok'
#endif

    return
  end subroutine init_swirljet

  subroutine postprocess_swirljet (ux1,uy1,uz1,pp3,phi1,ep1)

    use var, ONLY : nzmsize
    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

  end subroutine postprocess_swirljet

  subroutine visu_swirljet_init (visu_initialised)

    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D

    implicit none

    logical, intent(out) :: visu_initialised

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
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ur, utheta

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
