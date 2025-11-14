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

  subroutine boundary_conditions_swirljet (ux,uy,uz,phi)

    USE param
    USE variables

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux, uy, uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    call inflow (phi)
    call outflow (ux, uy, uz, phi)

    return
  end subroutine boundary_conditions_swirljet
  
  subroutine inflow (phi)

    USE param
    USE variables
    USE ibm_param

    implicit none

    integer  :: j, k, is
    real(mytype) :: y, z, r, theta
    real(mytype) :: um, alpha, delta
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    alpha = 100._mytype
    delta = theta_jet
    
    call random_number(bxo)
    call random_number(byo)
    call random_number(bzo)
    
    do k = 1, xsize(3)
       z = (k+xstart(3)-2)*dz-zlz/two
       do j = 1, xsize(2)
          if (istret.eq.0) y = (j+xstart(2)-2)*dy-yly/two
          if (istret.ne.0) y = yp(j+xstart(2)-1)-yly/two
          r = sqrt(y*y+z*z)
          theta = atan2(y,z)
          um = (one-erf((r-one)/delta))/two
          bxx1(j,k) = u1*(one-(alpha-one)/alpha/two*(one+erf((r-one)/delta))) &
                      +(bxo(j,k)-half)*um*inflow_noise
          bxy1(j,k) = u2*r/two*(one-erf((r-one)/delta))*cos(theta) &
                      +(byo(j,k)-half)*um*inflow_noise
          bxz1(j,k) = -u2*r/two*(one-erf((r-one)/delta))*sin(theta) &
                      +(bzo(j,k)-half)*um*inflow_noise
       enddo
    enddo

    if (iscalar.eq.1) then
       do is = 1, numscalar
          do k = 1,xsize(3)
             z = (k+xstart(3)-2)*dz-zlz/two
             do j = 1,xsize(2)
                if (istret.eq.0) y = (j+xstart(2)-2)*dy-yly/two
                if (istret.ne.0) y = yp(j+xstart(2)-1)-yly/two
                r = sqrt(y*y+z*z)
                phi(1,j,k,is) = cp(is)*(one-erf((r-one)/delta))/two
             enddo
          enddo
       enddo
    endif

    return
  end subroutine inflow

  subroutine outflow (ux, uy, uz, phi)

    USE param
    USE variables
    USE MPI
    USE ibm_param

    implicit none

    integer :: j,k,code
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx,uxmin,uxmax

    udx = one/dx; udy = one/dy; udz = one/dz;
    uddx = half/dx; uddy = half/dy; uddz = half/dz

    uxmax = -1609._mytype
    uxmin = 1609._mytype
    do k=1,xsize(3)
       do j=1,xsize(2)
          if (ux(nx-1,j,k).gt.uxmax) uxmax = ux(nx-1,j,k)
          if (ux(nx-1,j,k).lt.uxmin) uxmin = ux(nx-1,j,k)
       enddo
    enddo

    call MPI_ALLREDUCE(MPI_IN_PLACE,uxmax,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(MPI_IN_PLACE,uxmin,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)

    if (u1 == zero) then
       cx = (half*(uxmax+uxmin))*gdt(itr)*udx
    elseif (u1 == one) then
       cx = uxmax*gdt(itr)*udx
    elseif (u1 == two) then
       cx = u2*gdt(itr)*udx
    else
       cx = uxmax*gdt(itr)*udx
    endif

    do k = 1, xsize(3)
       do j = 1, xsize(2)
          bxxn(j,k) = ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
          bxyn(j,k) = uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
          bxzn(j,k) = uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
       enddo
    enddo

    if (iscalar==1) then
       if (u2==zero) then
          cx = (half*(uxmax+uxmin))*gdt(itr)*udx
       elseif (u2==one) then
          cx = uxmax*gdt(itr)*udx
       elseif (u2==two) then
          cx = u2*gdt(itr)*udx    !works better
       else
          !stop
          cx = uxmax*gdt(itr)*udx
       endif

       do k = 1, xsize(3)
          do j = 1, xsize(2)
             phi(nx,j,k,:) = phi(nx,j,k,:)-cx*(phi(nx,j,k,:)-phi(nx-1,j,k,:))
          enddo
       enddo
    endif

    if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) &
       write(*,*) "Outflow velocity ux nx=n min max=",real(uxmin,4),real(uxmax,4)

    return
  end subroutine outflow

  subroutine init_swirljet (ux1, uy1, uz1, phi1)

    USE decomp_2d_io
    USE variables
    USE param
    USE MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    integer :: i, j, k, ii, is, code
    real(mytype) :: y, z, r, theta
    real(mytype) :: um, alpha, delta

    ux1 = zero; uy1 = zero; uz1 = zero
    
    alpha = 100._mytype
    delta = theta_jet
    
    if (iin.ne.0) then
       call system_clock(count = code)
       if (iin.eq.2) code = 0
       call random_seed(size = ii)
       call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))

       call random_number(ux1)
       call random_number(uy1)
       call random_number(uz1)

       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                ux1(i,j,k) = init_noise*(ux1(i,j,k)-half)
                uy1(i,j,k) = init_noise*(uy1(i,j,k)-half)
                uz1(i,j,k) = init_noise*(uz1(i,j,k)-half)
             enddo
          enddo
       enddo

       !modulation of the random noise
       do k = 1, xsize(3)
          z = (k+xstart(3)-2)*dz-zlz/two
          do j = 1, xsize(2)
             if (istret.eq.0) y = (j+xstart(2)-2)*dy-yly/two
             if (istret.ne.0) y = yp(j+xstart(2)-1)-yly/two
             r = sqrt(y*y+z*z)
             um = (one-erf((r-one)/delta))/two
             do i = 1, xsize(1)
                ux1(i,j,k) = um*ux1(i,j,k)
                uy1(i,j,k) = um*uy1(i,j,k)
                uz1(i,j,k) = um*uz1(i,j,k)
             enddo
          enddo
       enddo
    endif

    !INIT FOR G AND U=MEAN FLOW + NOISE
    do k = 1, xsize(3)
       z = (k+xstart(3)-2)*dz-zlz/two
       do j = 1,xsize(2)
          if (istret.eq.0) y = (j+xstart(2)-2)*dy-yly/two
          if (istret.ne.0) y = yp(j+xstart(2)-1)-yly/two
          r = sqrt(y*y+z*z)
          theta = atan2(y,z)
          do i = 1,xsize(1)
             ux1(i,j,k) = ux1(i,j,k) + u1*(one-(alpha-one)/alpha/two*(one+erf((r-one)/delta)))
             uy1(i,j,k) = uy1(i,j,k) + u2*r/two*(one-erf((r-one)/delta))*cos(theta)
             uz1(i,j,k) = uz1(i,j,k) - u2*r/two*(one-erf((r-one)/delta))*sin(theta)
          enddo
       enddo
    enddo

    if (iscalar==1) then
       do k = 1, xsize(3)
          z = (k+xstart(3)-2)*dz-zlz/two
          do j = 1, xsize(2)
             if (istret.eq.0) y = (j+xstart(2)-2)*dy-yly/two
             if (istret.ne.0) y = yp(j+xstart(2)-1)-yly/two
             r = sqrt(y*y+z*z)
             do i = 1, xsize(1)
                do is = 1, numscalar
                   phi1(i,j,k,is) = cp(is)*(one-erf((r-one)/delta))/two
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
    !di1 = zero
    !di1(:,:,:)=sqrt(  (tf1(:,:,:)-th1(:,:,:))**2 &
    !                + (tg1(:,:,:)-tc1(:,:,:))**2 &
    !                + (tb1(:,:,:)-td1(:,:,:))**2)
    !call write_field(di1, ".", "vort", num, flush = .true.) ! Reusing temporary array, force flush

    !Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
    di1 = zero
    di1(:,:,:) = - half*(ta1(:,:,:)**2 + te1(:,:,:)**2 + ti1(:,:,:)**2) &
                 - td1(:,:,:) * tb1(:,:,:) &
                 - tg1(:,:,:) * tc1(:,:,:) &
                 - th1(:,:,:) * tf1(:,:,:)
    call write_field(di1, ".", "critq", num, flush = .true.) ! Reusing temporary array, force flush
    
    !lambda2
    !ta1=dudx te1=dudy ti1=dwdz
    !td1=0.5*(dudy+dvdx) tg1=0.5*(dudz+dwdz) th1=0.5*(dvdz+dwdy)
    !tb1=0.5*(dvdx-dudy) tc1=0.5*(dwdz-dudz) tf1=0.5*(dwdy-dvdz)
    di1 = zero
    td1(:,:,:) = half*(td1(:,:,:)+tb1(:,:,:)) !Sxy
    tg1(:,:,:) = half*(tg1(:,:,:)+tc1(:,:,:)) !Sxz
    th1(:,:,:) = half*(th1(:,:,:)+tf1(:,:,:)) !Syz
    tb1(:,:,:) = tb1(:,:,:) - td1(:,:,:) !Oxy
    tc1(:,:,:) = tc1(:,:,:) - tg1(:,:,:) !Oxz
    tf1(:,:,:) = tf1(:,:,:) - th1(:,:,:) !Oyz
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             m11 = ta1(i,j,k)*ta1(i,j,k)+td1(i,j,k)*td1(i,j,k)+tg1(i,j,k)*tg1(i,j,k) &
                  -tb1(i,j,k)*tb1(i,j,k)-tc1(i,j,k)*tc1(i,j,k)
             m12 = ta1(i,j,k)*td1(i,j,k)+td1(i,j,k)*te1(i,j,k)+tg1(i,j,k)*th1(i,j,k) &
                  -tc1(i,j,k)*tf1(i,j,k)
             m13 = ta1(i,j,k)*tg1(i,j,k)+td1(i,j,k)*th1(i,j,k)+tg1(i,j,k)*ti1(i,j,k) &
                  +tb1(i,j,k)*tf1(i,j,k)
             m22 = td1(i,j,k)*td1(i,j,k)+te1(i,j,k)*te1(i,j,k)+th1(i,j,k)*th1(i,j,k) &
                  -tb1(i,j,k)*tb1(i,j,k)-tf1(i,j,k)*tf1(i,j,k)
             m23 = td1(i,j,k)*tg1(i,j,k)+te1(i,j,k)*th1(i,j,k)+th1(i,j,k)*ti1(i,j,k) &
                  -tb1(i,j,k)*tc1(i,j,k)
             m33 = tg1(i,j,k)*tg1(i,j,k)+th1(i,j,k)*th1(i,j,k)+ti1(i,j,k)*ti1(i,j,k) &
                  -tc1(i,j,k)*tc1(i,j,k)-tf1(i,j,k)*tf1(i,j,k)
             trace = m11+m22+m33
             det = m11*(m22*m33-m23*m23)-m12*(m12*m33-m23*m13)+m13*(m12*m23-m22*m13)
             b = m11*m22+m11*m33+m22*m33-m12*m12-m13*m13-m23*m23
             q = trace*trace/nine-b/three
             r = trace*trace*trace/twentyseven-trace*b/six-det/two
             
             eig1 = trace/three; eig2 = eig1; eig3 = eig1;
             if (q<0.0) then
                tmp = sqrt(-q*q*q);
                if (abs(r) <= tmp) then
                   theta = acos(r/tmp)
                   eig1 = trace/three+two*sqrt(-q)*cos(theta/three)
                   eig2 = trace/three+two*sqrt(-q)*cos((theta+twopi)/three)
                   eig3 = trace/three+two*sqrt(-q)*cos((theta-twopi)/three)
                   if (eig1>eig2) call swap(eig1,eig2)
                   if (eig2>eig3) call swap(eig2,eig3)
                   if (eig1>eig2) call swap(eig1,eig2)
                endif
             endif
             
             di1(i,j,k) = eig2
          enddo
       enddo
    enddo
    call write_field(di1, ".", "lambda2", num, flush = .true.) ! Reusing temporary array, force flush
    
    do k = 1, xsize(3)
       z = (k+xstart(3)-2)*dz-zlz/two
       do j=1,xsize(2)
          if (istret.eq.0) y = (j+xstart(2)-2)*dy-yly/two
          if (istret.ne.0) y = yp(j+xstart(2)-1)-yly/two
          theta = atan2(y,z)
          do i = 1, xsize(1)
             ur(i,j,k) = uz1(i,j,k)*cos(theta) + uy1(i,j,k)*sin(theta)
             utheta(i,j,k) = -uz1(i,j,k)*sin(theta) + uy1(i,j,k)*cos(theta)
          enddo
       enddo
    enddo
    call write_field(ur, ".", "ur", num)
    call write_field(utheta, ".", "utheta", num)

  end subroutine visu_swirljet
  
  subroutine swap(a,b)
  
    implicit none
    real(mytype), intent(inout) :: a,b
    real(mytype) :: Temp

    Temp = a
    a = b
    b = Temp

  end subroutine swap

end module swirljet
