!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module rbc

  USE decomp_2d_constants
  USE decomp_2d_mpi
  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  PRIVATE !! All functions/subroutines private by default
  PUBLIC :: geomcomplex_rbc, init_rbc, boundary_conditions_rbc, postprocess_rbc

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

  subroutine init_rbc (ux1,uy1,uz1,ep1,phi1)

    use decomp_2d_io
    use variables
    use param
    use ibm_param
    use MPI

    implicit none

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    integer :: i, j, k, is, ii, code
    real(mytype) :: x, y, z, r, um, disturb

    ux1=zero; uy1=zero; uz1=zero

    if (iin.ne.0) then
       call system_clock(count = code)
       if (iin.eq.2) code = 0
       call random_seed(size = ii)
       call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /))
    endif

    if (iscalar.eq.1) then
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
                                phi1(i,j,k,is) = one - y/yly + um*disturb*init_noise
                            else
                                phi1(i,j,k,is) = zero
                            endif
                        else
                            phi1(i,j,k,is) = one - y/yly + um*disturb*init_noise
                        endif
                    enddo
                enddo
            enddo
        enddo
    endif

    if (nrank  ==  0) write(*,*) '# init end ok'

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

    return
  end subroutine boundary_conditions_rbc
  !########################################################################

  subroutine postprocess_rbc(ux1,uy1,uz1,pp3,phi1,ep1)

    use var, ONLY : nzmsize

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

  end subroutine postprocess_rbc
  !########################################################################

end module rbc
