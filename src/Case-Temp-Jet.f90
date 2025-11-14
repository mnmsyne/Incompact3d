!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module tempjet

  USE decomp_2d_constants
  USE decomp_2d_mpi
  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  PRIVATE !! All functions/subroutines private by default
  PUBLIC :: init_tempjet, postprocess_tempjet, visu_tempjet, visu_tempjet_init

  ! M. Hayashi, T. Watanabe, and K. Nagata,
  ! Characteristics of small-scale shear layers in a temporally evolving turbulent planar jet,
  ! J. Fluid Mech. 920, A38 (2021).

contains

  subroutine init_tempjet (ux1, uy1, uz1, ep1, phi1)

    use decomp_2d_io
    use param
    use variables
    use MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    integer :: i, j, k, ii, is, code
    real(mytype) :: x, y, z, r, um

    integer :: Nxi, Nyi, Nzi
    real(mytype) :: Delta_in
    real(mytype), allocatable, dimension(:,:,:) :: ux_coarse, uy_coarse, uz_coarse
    real(mytype), allocatable, dimension(:) :: x_coarse, y_coarse, z_coarse

   !  zeromach=one
   !  do while ((one + zeromach / two) .gt. one)
   !     zeromach = zeromach/two
   !  end do
   !  zeromach = ten*zeromach
    
    ux1 = zero; uy1 = zero; uz1 = zero

    Delta_in = 0.02_mytype

    if (iin.eq.1) then !generation of a random noise     
       call system_clock(count = code)
       call random_seed(size = ii)
       call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))
       
       call random_number(ux1)
       call random_number(uy1)
       call random_number(uz1)
       
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                ux1(i,j,k) = init_noise*(ux1(i,j,k)-0.5)
                uy1(i,j,k) = init_noise*(uy1(i,j,k)-0.5)
                uz1(i,j,k) = init_noise*(uz1(i,j,k)-0.5)
             enddo
          enddo
       enddo
       
       !modulation of the random noise
       do k = 1, xsize(3)
          z = (k+xstart(3)-2)*dz-zlz/two
          do j = 1, xsize(2)
             if (istret.eq.0) y = (j+xstart(2)-2)*dy-yly/two
             if (istret.ne.0) y = yp(j+xstart(2)-1)-yly/two
             um = half+half*tanh((one-two*sqrt(y*y))/four/theta_jet)
             do i = 1, xsize(1)
                ux1(i,j,k) = um*ux1(i,j,k)
                uy1(i,j,k) = um*uy1(i,j,k)
                uz1(i,j,k) = um*uz1(i,j,k)
             enddo
          enddo
       enddo
    endif

    if (iin.eq.2) then
       Nxi = ceiling(xlx / Delta_in) + 1
       Nyi = ceiling(yly / Delta_in) + 1
       Nzi = ceiling(zlz / Delta_in) + 1

       allocate(x_coarse(Nxi), y_coarse(Nyi), z_coarse(Nzi))
       allocate(ux_coarse(Nxi, Nyi, Nzi))
       allocate(uy_coarse(Nxi, Nyi, Nzi))
       allocate(uz_coarse(Nxi, Nyi, Nzi))

       do i = 1, Nxi
          x_coarse(i) = (i-1) * Delta_in
       end do
       do j = 1, Nyi
          y_coarse(j) = (j-1) * Delta_in - yly / 2.
       end do
       do k = 1, Nzi
          z_coarse(k) = (k-1) * Delta_in
       end do

       call system_clock(count=code)
       call random_seed(size = ii)
       call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))
      
       call random_number(ux_coarse)
       call random_number(uy_coarse)
       call random_number(uz_coarse)

       ux_coarse = init_noise * (ux_coarse - 0.5_mytype)
       uy_coarse = init_noise * (uy_coarse - 0.5_mytype)
       uz_coarse = init_noise * (uz_coarse - 0.5_mytype)

       do k = 1, xsize(3)
         z = (k + xstart(3) - 2) * dz
         do j = 1, xsize(2)
            if (istret.eq.0) y = (j + xstart(2) - 2) * dy -yly / 2.
            if (istret.ne.0) y = yp(j + xstart(2) - 1) - yly / 2.
            do i = 1, xsize(1)
                x = (i + xstart(1) - 2) * dx
                ux1(i,j,k) = interp3(x_coarse, y_coarse, z_coarse, ux_coarse, x, y, z, Nxi, Nyi, Nzi)
                uy1(i,j,k) = interp3(x_coarse, y_coarse, z_coarse, uy_coarse, x, y, z, Nxi, Nyi, Nzi)
                uz1(i,j,k) = interp3(x_coarse, y_coarse, z_coarse, uz_coarse, x, y, z, Nxi, Nyi, Nzi)
             enddo
          enddo
       enddo

       deallocate(x_coarse, y_coarse, z_coarse)
       deallocate(ux_coarse, uy_coarse, uz_coarse)
       
       !modulation of the random noise
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             if (istret.eq.0) y = (j+xstart(2)-2)*dy-yly/two
             if (istret.ne.0) y = yp(j+xstart(2)-1)-yly/two
             um = half + half * tanh((one - two * sqrt(y * y)) / four / theta_jet)
             do i = 1 ,xsize(1)
                ux1(i,j,k) = um*ux1(i,j,k)
                uy1(i,j,k) = um*uy1(i,j,k)
                uz1(i,j,k) = um*uz1(i,j,k)
             enddo
          enddo
       enddo
    endif
       
    do k = 1, xsize(3)
       z = (k+xstart(3)-2)*dz-zlz/two
       do j = 1, xsize(2)
          if (istret.eq.0) y = (j+xstart(2)-2)*dy-yly/two
          if (istret.ne.0) y = yp(j+xstart(2)-1)-yly/two
          do i = 1, xsize(1)
             !! Set mean field
             ux1(i, j, k) = ux1(i, j, k) + half * (u1 + u2) &
                  + half * (u1 - u2) * tanh((one - two * sqrt(y ** two)) / four / theta_jet)
             uy1(i, j, k) = uy1(i, j, k)
             uz1(i, j, k) = uz1(i, j, k)
          enddo
       enddo
    enddo
    
    if (iscalar.eq.1) then
       do k = 1, xsize(3)
          z = (k+xstart(3)-2)*dz-zlz/two
          do j = 1, xsize(2)
             if (istret.eq.0) y = (j+xstart(2)-2)*dy-yly/two
             if (istret.ne.0) y = yp(j+xstart(2)-1)-yly/two
             do i = 1, xsize(1)
                do is = 1, numscalar
                   phi1(i,j,k,is) = cp(is) * (half &
                       + half * tanh((one - two * sqrt(y ** two)) / four / theta_jet))
                enddo
             enddo
          enddo
       enddo
    endif

#ifdef DEBG
    if (nrank  ==  0) write(*,*) '# init end ok'
#endif

    return

    contains

    function interp3(x, y, z, f, xi, yi, zi, nx, ny, nz) result(f_interp)
      integer, intent(in) :: nx, ny, nz
      real(mytype), intent(in) :: x(nx), y(ny), z(nz)
      real(mytype), intent(in) :: f(nx, ny, nz)
      real(mytype), intent(in) :: xi, yi, zi
      real(mytype) :: f_interp
      
      integer :: i, j, k
      real(mytype) :: xd, yd, zd
      real(mytype) :: c00, c01, c10, c11, c0, c1
      
      i = 1
      do while (i < nx .and. x(i+1) < xi)
         i = i + 1
      end do
      
      j = 1
      do while (j < ny .and. y(j+1) < yi)
         j = j + 1
      end do
      
      k = 1
      do while (k < nz .and. z(k+1) < zi)
         k = k + 1
      end do
      
      if (i < nx .and. j < ny .and. k < nz) then
         xd = (xi - x(i)) / (x(i+1) - x(i))
         yd = (yi - y(j)) / (y(j+1) - y(j))
         zd = (zi - z(k)) / (z(k+1) - z(k))
         
         c00 = f(i,j,k) * (1-xd) + f(i+1,j,k) * xd
         c01 = f(i,j,k+1) * (1-xd) + f(i+1,j,k+1) * xd
         c10 = f(i,j+1,k) * (1-xd) + f(i+1,j+1,k) * xd
         c11 = f(i,j+1,k+1) * (1-xd) + f(i+1,j+1,k+1) * xd
         
         c0 = c00 * (1-yd) + c10 * yd
         c1 = c01 * (1-yd) + c11 * yd
         
         f_interp = c0 * (1-zd) + c1 * zd
      else
         f_interp = 0.0_mytype
      endif
      
    end function interp3

  end subroutine init_tempjet
  !########################################################################

  subroutine postprocess_tempjet(ux1, uy1, uz1, pp3, phi1, ep1)

    use var, ONLY : nzmsize
    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

  end subroutine postprocess_tempjet
  !########################################################################

  subroutine visu_tempjet_init (visu_initialised)

    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "critq", 1, 0, output2D, mytype)

    visu_initialised = .true.
    
  end subroutine visu_tempjet_init
  !########################################################################
  
  subroutine visu_tempjet(ux1, uy1, uz1, pp3, phi1, ep1, num)

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

  end subroutine visu_tempjet
  !########################################################################
  
end module tempjet
