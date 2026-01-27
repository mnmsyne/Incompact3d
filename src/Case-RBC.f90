!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module rbc

  USE decomp_2d_constants
  USE decomp_2d_mpi
  USE decomp_2d
  USE variables
  USE param

  implicit none

  PRIVATE !! All functions/subroutines private by default
  PUBLIC :: geomcomplex_rbc, &
            init_rbc, &
            boundary_conditions_rbc, &
            postprocess_rbc

contains

  subroutine geomcomplex_rbc (epsi, nxi, nxf, ny, nyi, nyf, nzi, nzf, dx, yp, dz, remp)

    use MPI
    use param, only : zero, one, two, xlx, zlz
    use ibm_param

    implicit none

    integer :: nxi, nxf, ny, nyi, nyf, nzi, nzf
    real(mytype), dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
    real(mytype), dimension(ny) :: yp
    real(mytype) :: dx, dz, remp, tol

    integer :: i, j, k, code, ierror
    real(mytype) :: r, xm, zm
    
    epsi(:,:,:) = zero
    tol = 1e-15

    !safety check
    if (nrank == 0) then
       write(*,*) 'WARNING: IN PROGRESS - cylindrical RBC'
       if (ra.le.0) then
           write(*,*) 'SIMULATION IS STOPPED!'
           write(*,*) 'Please specify a valid value for the radius in input.i3d (ra)'
           call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
       endif
    endif

    do k = nzi, nzf
        zm = real(k-1,mytype)*dz-zlz/two
        do j = nyi, nyf
            do i = nxi, nxf
                xm = real(i-1,mytype)*dx-xlx/two
                r = sqrt(xm*xm+zm*zm)
                if (r.gt.ra) then
                   epsi(i,j,k) = remp
                elseif (abs(r-ra).lt.tol) then
                    epsi(i,j,k) = remp
                endif
            enddo
        enddo
    enddo

    return
  end subroutine geomcomplex_rbc
  !########################################################################

  subroutine init_rbc (ux1, uy1, uz1, ep1, phi1)

    use decomp_2d_io
    use variables
    use param
    use ibm_param
    use MPI

    implicit none

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    integer :: i, j, k, is, ii, code, ierror
    real(mytype) :: x, y, z, r, um, disturb

    if (iscalar.eq.0) then
        if (nrank.eq.0) then
            write(*,*) 'SIMULATION IS STOPPED!'
            write(*,*) 'Please specify a scalar for RBC'
            call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
        endif
    endif

    ux1=zero; uy1=zero; uz1=zero

    if (iin.ne.0) then
       call system_clock(count = code)
       if (iin.eq.2) code = 0
       call random_seed(size = ii)
       call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))
    endif

    do is = 1, numscalar
        do k = 1, xsize(3)
            z = real(k+xstart(3)-2,mytype)*dz-zlz/two
            do j = 1, xsize(2)
                if (istret.eq.0) y = real(j+xstart(2)-2,mytype)*dy
                if (istret.ne.0) y = yp(j+xstart(2)-1)
                um = sin(PI*y/yly)
                do i = 1, xsize(1)
                    x = real(i+xstart(1)-2,mytype)*dx-xlx/two
                    r = sqrt(x*x+z*z)

                    if (iin.eq.0) then
                        disturb = cos(two*PI*x/xlx)/sixteen + cos(four*PI*x/xlx)/eight + cos(eight*PI*x/xlx) &
                                + cos(two*PI*z/zlz)/sixteen + cos(four*PI*z/zlz)/eight + cos(eight*PI*z/zlz)
                    else
                        call random_number(disturb)
                        disturb = two*disturb - one
                    endif

                    if (iibm.ne.0) then
                        if (r.le.ra .and. ep1(i,j,k).eq.0) then
                            phi1(i,j,k,is) = one - y/yly + init_noise*disturb*um*cos(half*PI*r/ra)
                        else
                            phi1(i,j,k,is) = zero
                        endif
                    else
                        phi1(i,j,k,is) = one - y/yly + init_noise*disturb*um
                    endif
                enddo
            enddo
        enddo
    enddo

    if (nrank.eq.0) write(*,*) '# init end ok'

    return
  end subroutine init_rbc
  !########################################################################

  subroutine boundary_conditions_rbc (ux, uy, uz, phi)

    use param
    use variables
    use MPI

    implicit none

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux, uy, uz
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    integer :: is

    if (iscalar.ne.0) then
        ! y = 0
        if (nclyS1==2 .and. xstart(2)==1) then
            do is = 1, numscalar
                phi(:,1,:,is) = one
            enddo
        endif
        ! y = Ly
        if (nclySn==2 .and. xend(2)==ny) then
            do is = 1, numscalar
                phi(:,xsize(2),:,is) = zero
            enddo
        endif
    endif

  end subroutine boundary_conditions_rbc
  !########################################################################

  subroutine postprocess_rbc(ux1, uy1, uz1, pp3, phi1, ep1)

    use param, only : zero, itime, iprocessing, t, npress
    use variables, only : numscalar
    use var, only : nzmsize

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

    integer :: FS
    character(len=100) :: fileformat
    character(len=1), parameter :: NL=char(10) !new line character
    logical, save :: init = .FALSE.

    integer :: i, j, k, is
    real(mytype) :: nu_bot, nu_top, nu_vol

    if (.not.init) then
        if (nrank == 0) write(*,*) '# postprocess RBC initialized'
        init = .TRUE.  
    endif

    if (mod(itime,iprocessing).ne.0) return

    nu_bot = zero; nu_top = zero; nu_vol = zero

    call wall_nusselt_rbc(phi1, nu_bot, nu_top)
    call volume_nusselt_rbc(uy1, phi1, nu_vol)

    if (nrank .eq. 0) then
        FS = 1+3 !Number of columns
        write(fileformat, '( "(",I4,"(E14.6),A)" )' ) FS
        FS = FS*14+1 !Line width
        open(67, file='./out/nusselt.txt', status='unknown', form='formatted', &
            access='direct', recl=FS)
        write(67, fileformat, rec=itime/iprocessing+1) &
            t, &             !time
            nu_bot, &        !nusselt number at bottom wall
            nu_top, &        !nusselt number at top wall
            nu_vol, &        !volume nusselt number
            NL
        close(67)
    end if

  end subroutine postprocess_rbc
  !########################################################################

  subroutine wall_nusselt_rbc(phi1, nu_bot1, nu_top1)

    USE decomp_2d_io
    USE MPI

    use variables, only : nx, ny, nz, yp, numscalar
    use param, only : istret, dy
    use param, only : zero, two, six, nine, eleven, eighteen
    use var, only : phi2

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(out) :: nu_bot1, nu_top1

    real(mytype), dimension(ysize(1), ysize(3)) :: wallnu2
    integer :: i, k, is, code
    real(mytype) :: nu_bot, nu_top
    real(mytype) :: dy0, dy1

    nu_bot = zero; nu_top = zero
    nu_bot1 = zero; nu_top1 = zero

    is = 1
    call transpose_x_to_y(phi1(:,:,:,is), phi2(:,:,:,is))

    !bottom wall
    wallnu2 = zero
    do k = 1, ysize(3)
        do i = 1, ysize(1)
            if (istret.eq.0) then
                dy0 = dy
                wallnu2(i,k) = -(- eleven * phi2(i,1,k,is) &
                                + eighteen * phi2(i,2,k,is) &
                                - nine * phi2(i,3,k,is) &
                                + two * phi2(i,4,k,is)) / six / dy0
            else
                dy0 = yp(2) - yp(1)
                dy1 = yp(3) - yp(2)
                wallnu2(i,k) = -(- phi2(i,1,k,is) * (two * dy0 + dy1) / dy0 / (dy0 + dy1) &
                                + phi2(i,2,k,is) * (dy0 + dy1) / dy0 / dy1 &
                                - phi2(i,3,k,is) * dy0 / dy1 / (dy0 + dy1))
            endif
            nu_bot = nu_bot + wallnu2(i,k)
        enddo
    enddo

    !top wall
    wallnu2 = zero
    do k = 1, ysize(3)
        do i = 1, ysize(1)
            if (istret.eq.0) then
                dy0 = dy
                wallnu2(i,k) = -(eleven * phi2(i,ny,k,is) &
                                - eighteen * phi2(i,ny-1,k,is) &
                                + nine * phi2(i,ny-2,k,is) &
                                - two * phi2(i,ny-3,k,is)) / six / dy0
            else
                dy0 = yp(ny) - yp(ny-1)
                dy1 = yp(ny-1) - yp(ny-2)
                wallnu2(i,k) = -(phi2(i,ny,k,is) * (two * dy0 + dy1) / dy0 / (dy0 + dy1) &
                                - phi2(i,ny-1,k,is) * (dy0 + dy1) / dy0 / dy1 &
                                + phi2(i,ny-2,k,is) * dy0 / dy1 / (dy0 + dy1))
            endif
            nu_top = nu_top + wallnu2(i,k)
        enddo
    enddo

    call MPI_REDUCE(nu_bot,nu_bot1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
    call MPI_REDUCE(nu_top,nu_top1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

    if (nrank .eq. 0) then
        nu_bot1 = nu_bot1 / real(nx*nz,mytype)
        nu_top1 = nu_top1 / real(nx*nz,mytype)
    endif

  end subroutine wall_nusselt_rbc
  !########################################################################

  subroutine volume_nusselt_rbc(uy1, phi1, nu_vol1)

    USE decomp_2d_io
    USE MPI

    use param, only : u1, u2, istret, dx, dy, dz, yly
    use param, only : zero, one, half
    use variables, only : derys
    use variables, only : nx, ny, nz, yp, numscalar
    use var, only : phi2, uy2
    ! use var, only : di2, sy, ffys, fsys, fwys, ppy

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: uy1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(out) :: nu_vol1

    ! real(mytype), dimension(ysize(1),ysize(2),ysize(3)) :: dphiy2
    integer :: i, j, k, is, code
    real(mytype) :: nu_vol, theta, weight, y

    ! rayleigh = u1; prandtl = u2
    nu_vol = zero; nu_vol1 = zero

    is = 1
    call transpose_x_to_y(uy1, uy2)
    call transpose_x_to_y(phi1(:,:,:,is), phi2(:,:,:,is))
    ! call deryS (dphiy2, phi2(:,:,:,is), di2, sy, ffyS, fsyS, fwyS, ppy, ysize(1), ysize(2), ysize(3), 1, zero)

    if (istret.eq.0) then
        do k = 1, ysize(3)
            do j = 1, ysize(2)
                y = real(j+ystart(2)-2,mytype)*dy
                do i = 1, ysize(1)
                    theta = phi2(i,j,k,is) - (one - y/yly)
                    nu_vol = nu_vol + uy2(i,j,k) * theta
                enddo
            enddo
        enddo

        call MPI_REDUCE(nu_vol,nu_vol1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

        if (nrank .eq. 0) then
            nu_vol1 = nu_vol1 / (real(nx*(ny-1)*nz,mytype))
            nu_vol1 = one + sqrt(u1*u2) * nu_vol1
        endif
    else
        do k = 1, ysize(3)
            do j = 1, ysize(2)
                y = yp(j+ystart(2)-1)
                if (j == 1) then
                    weight = half * (yp(2) - yp(1))
                elseif (j == ny) then
                    weight = half * (yp(ny) - yp(ny-1))
                else
                    weight = half * (yp(j+1) - yp(j-1))
                endif
                do i = 1, ysize(1)
                    theta = phi2(i,j,k,is) - (one - y/yly)
                    nu_vol = nu_vol + weight * uy2(i,j,k) * theta
                enddo
            enddo
        enddo

        call MPI_REDUCE(nu_vol,nu_vol1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

        if (nrank .eq. 0) then
            nu_vol1 = nu_vol1 / yly / (real(nx*nz,mytype))
            nu_vol1 = one + sqrt(u1*u2) * nu_vol1
        endif
    endif

  end subroutine volume_nusselt_rbc
  !########################################################################

end module rbc
