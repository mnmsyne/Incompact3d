!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module rectosc

  USE decomp_2d_constants
  USE decomp_2d_mpi
  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: init_rectosc, boundary_conditions_rectosc, postprocess_rectosc, &
            geomcomplex_rectosc

contains

  subroutine rectosc_motion(xc, yc)

    use param, only : zero, one, two
    use ibm_param

    implicit none

    real(mytype), intent(out) :: xc, yc
    real(mytype) :: amp, tau, pi_local

    pi_local = acos(-one)
    tau = zero
    if (t /= zero) tau = t - real(ifirst, mytype) * dt

    amp = rect_amp
    if (amp <= zero .and. omega > zero) amp = pi_local / (two * omega)

    xc = cex
    yc = cey
    ubcx = zero
    ubcy = zero
    ubcz = zero

    if (omega > zero) then
       xc = cex + amp * sin(omega * tau)
       ubcx = amp * omega * cos(omega * tau)
    endif

  end subroutine rectosc_motion

  subroutine geomcomplex_rectosc(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,dx,yp,remp)

    use param, only : zero, half
    use ibm_param

    implicit none

    integer :: nxi,nxf,ny,nyi,nyf,nzi,nzf
    real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
    real(mytype),dimension(ny) :: yp
    real(mytype) :: dx, remp
    integer :: i,j,k
    real(mytype) :: xm,ym,xc,yc,xhalf,yhalf

    call rectosc_motion(xc, yc)

    epsi(:,:,:) = zero
    xhalf = half * rect_length
    yhalf = half * rect_width

    do k=nzi,nzf
       do j=nyi,nyf
          ym = yp(j)
          do i=nxi,nxf
             xm = real(i-1,mytype) * dx
             if (abs(xm - xc) <= xhalf .and. abs(ym - yc) <= yhalf) then
                epsi(i,j,k) = remp
             endif
          enddo
       enddo
    enddo

  end subroutine geomcomplex_rectosc

  subroutine init_rectosc(ux1,uy1,uz1,phi1)

    use cyl, only : init_cyl

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    call init_cyl(ux1, uy1, uz1, phi1)

  end subroutine init_rectosc

  subroutine boundary_conditions_rectosc(ux,uy,uz,phi)

    use cyl, only : boundary_conditions_cyl

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    call boundary_conditions_cyl(ux, uy, uz, phi)

  end subroutine boundary_conditions_rectosc

  subroutine postprocess_rectosc(ux1,uy1,uz1,ep1)

    implicit none

    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1

  end subroutine postprocess_rectosc

end module rectosc
