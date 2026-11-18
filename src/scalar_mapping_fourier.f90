!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module scalar_mapping_fourier

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : decomp_2d_abort,real_type
  use decomp_2d, only : xsize,ysize,zsize,xstart, &
       transpose_x_to_y,transpose_y_to_z,transpose_z_to_y,transpose_y_to_x
  use variables, only : nx,ny,nz
  use param, only : dx,dy,dz
  use scalar_flux, only : scalar_cell_width
  use mpi
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
#ifdef X3D_SCALAR_MAPPING_MKL
  use MKL_DFTI
#endif
  implicit none
  private
  public :: scalar_fourier_available,scalar_fourier_solve,scalar_fourier_finalize

  logical, save :: initialized=.false.
  complex(mytype), allocatable, save :: sx(:,:,:),sy(:,:,:),sz(:,:,:)
  real(mytype), allocatable, save :: inverse_pivot(:,:,:),width(:)
#ifdef X3D_SCALAR_MAPPING_MKL
  type(DFTI_DESCRIPTOR), pointer, save :: plan_y=>null(),plan_z=>null()
#endif

contains

  logical function scalar_fourier_available()
    scalar_fourier_available=.false.
#ifdef X3D_SCALAR_MAPPING_MKL
    scalar_fourier_available=.true.
#endif
  end function scalar_fourier_available

  subroutine initialize_fourier()
    integer :: i,j,k,m,n,count
    real(mytype) :: mu,pivot,pi,previous_upper

    if (initialized) return
    if (.not.scalar_fourier_available()) &
         call decomp_2d_abort(1,'scalar_mapping_solver=fourier requires the MKL mapping backend')
    allocate(sx(xsize(1),xsize(2),xsize(3)),sy(ysize(1),ysize(2),ysize(3)), &
         sz(zsize(1),zsize(2),zsize(3)))
    allocate(inverse_pivot(nx-1,xsize(2),xsize(3)),width(nx-1))
    do i=1,nx-1
       width(i)=scalar_cell_width(i+1)
    enddo
    pi=acos(-1._mytype)
    ! W*Lx has off-diagonals -1/dx, including the unequal end volumes.
    ! Factor each shifted real tridiagonal matrix once; complex RHS are solved
    ! in place. No rank replicates the physical or spectral global field.
    do k=1,xsize(3)
       n=xstart(3)+k-2
       do j=1,xsize(2)
          m=xstart(2)+j-2
          mu=4._mytype*sin(pi*real(m,mytype)/real(ny,mytype))**2/dy**2 &
               +4._mytype*sin(pi*real(n,mytype)/real(nz,mytype))**2/dz**2
          do i=1,nx-1
             count=2
             if (i==1.or.i==nx-1) count=1
             pivot=real(count,mytype)/dx+mu*width(i)
             if (m==0.and.n==0.and.i==1) then
                ! Pin only lambda at the first ACTIVE x node. The rest of the
                ! transverse zero mode remains a nontrivial 1D Neumann solve.
                pivot=1._mytype
             elseif (i>1) then
                previous_upper=-1._mytype/dx
                if (m==0.and.n==0.and.i==2) previous_upper=0._mytype
                pivot=pivot+previous_upper*inverse_pivot(i-1,j,k)/dx
             endif
             if (.not.ieee_is_finite(pivot).or.pivot<=tiny(pivot)) &
                  call decomp_2d_abort(1,'Invalid scalar Fourier mapping tridiagonal pivot')
             inverse_pivot(i,j,k)=1._mytype/pivot
          enddo
       enddo
    enddo
#ifdef X3D_SCALAR_MAPPING_MKL
    call make_plan(plan_y,ny,ysize(1),ysize(1))
    call make_plan(plan_z,nz,zsize(1)*zsize(2),zsize(1)*zsize(2))
#endif
    initialized=.true.
  end subroutine initialize_fourier

  subroutine scalar_fourier_solve(rhs,potential)
    real(mytype), intent(in) :: rhs(xsize(1),xsize(2),xsize(3))
    real(mytype), intent(out) :: potential(xsize(1),xsize(2),xsize(3))
    real(mytype) :: local(2),global(2)
    integer :: i,j,k,m,n,ierr,bad,global_bad

    call initialize_fourier()
    sx=cmplx(rhs,0._mytype,kind=mytype)
    sx(1,:,:)=cmplx(0._mytype,0._mytype,kind=mytype)
    call transpose_x_to_y(sx,sy)
#ifdef X3D_SCALAR_MAPPING_MKL
    do k=1,ysize(3)
       call transform(plan_y,sy(:,:,k),.true.)
    enddo
#endif
    call transpose_y_to_z(sy,sz)
#ifdef X3D_SCALAR_MAPPING_MKL
    call transform(plan_z,sz,.true.)
#endif
    call transpose_z_to_y(sz,sy)
    call transpose_y_to_x(sy,sx)
    do k=1,xsize(3)
       n=xstart(3)+k-2
       do j=1,xsize(2)
          m=xstart(2)+j-2
          do i=1,nx-1
             sx(i+1,j,k)=width(i)*sx(i+1,j,k)
          enddo
          if (m==0.and.n==0) sx(2,j,k)=cmplx(0._mytype,0._mytype,kind=mytype)
          do i=2,nx-1
             sx(i+1,j,k)=sx(i+1,j,k)+sx(i,j,k)*inverse_pivot(i-1,j,k)/dx
          enddo
          sx(nx,j,k)=sx(nx,j,k)*inverse_pivot(nx-1,j,k)
          do i=nx-2,1,-1
             if (m==0.and.n==0.and.i==1) then
                sx(2,j,k)=cmplx(0._mytype,0._mytype,kind=mytype)
             else
                sx(i+1,j,k)=(sx(i+1,j,k)+sx(i+2,j,k)/dx)*inverse_pivot(i,j,k)
             endif
          enddo
       enddo
    enddo
    call transpose_x_to_y(sx,sy)
    call transpose_y_to_z(sy,sz)
#ifdef X3D_SCALAR_MAPPING_MKL
    call transform(plan_z,sz,.false.)
#endif
    call transpose_z_to_y(sz,sy)
#ifdef X3D_SCALAR_MAPPING_MKL
    do k=1,ysize(3)
       call transform(plan_y,sy(:,:,k),.false.)
    enddo
#endif
    call transpose_y_to_x(sy,sx)
    potential=real(sx,mytype)
    potential(1,:,:)=0._mytype
    bad=0
    if (.not.all(ieee_is_finite(potential)).or..not.all(ieee_is_finite(aimag(sx)))) bad=1
    call MPI_ALLREDUCE(bad,global_bad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (ierr/=MPI_SUCCESS) call decomp_2d_abort(1,'Scalar Fourier mapping finite-check reduction failed')
    if (global_bad/=0) call decomp_2d_abort(1,'Nonfinite scalar Fourier mapping solution')
    local=[maxval(abs(aimag(sx(2:nx,:,:)))),maxval(abs(potential(2:nx,:,:)))]
    call MPI_ALLREDUCE(local,global,2,real_type,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (ierr/=MPI_SUCCESS) call decomp_2d_abort(1,'Scalar Fourier mapping norm reduction failed')
    if (.not.all(ieee_is_finite(global)).or.global(1)>8192._mytype*epsilon(1._mytype)*max(1._mytype,global(2))) &
         call decomp_2d_abort(1,'Nonreal or nonfinite scalar Fourier mapping solution')
    ! The caller fixes the volume-weighted gauge and checks the original
    ! physical-space operator, not just the transformed equations.
  end subroutine scalar_fourier_solve

#ifdef X3D_SCALAR_MAPPING_MKL
  subroutine check_fft(status)
    integer, intent(in) :: status
    if (status/=0) call decomp_2d_abort(1,'MKL error in independent scalar mapping FFT plan')
  end subroutine check_fft

  subroutine make_plan(plan,length,batch,stride)
    type(DFTI_DESCRIPTOR), pointer, intent(inout) :: plan
    integer, intent(in) :: length,batch,stride
    integer :: precision
    precision=DFTI_SINGLE
    if (mytype==kind(1.d0)) precision=DFTI_DOUBLE
    call check_fft(DftiCreateDescriptor(plan,precision,DFTI_COMPLEX,1,length))
    call check_fft(DftiSetValue(plan,DFTI_PLACEMENT,DFTI_INPLACE))
    call check_fft(DftiSetValue(plan,DFTI_NUMBER_OF_TRANSFORMS,batch))
    call check_fft(DftiSetValue(plan,DFTI_INPUT_DISTANCE,1))
    call check_fft(DftiSetValue(plan,DFTI_OUTPUT_DISTANCE,1))
    call check_fft(DftiSetValue(plan,DFTI_INPUT_STRIDES,[0,stride]))
    call check_fft(DftiSetValue(plan,DFTI_OUTPUT_STRIDES,[0,stride]))
    call check_fft(DftiSetValue(plan,DFTI_BACKWARD_SCALE,1._mytype/real(length,mytype)))
    call check_fft(DftiCommitDescriptor(plan))
  end subroutine make_plan

  subroutine transform(plan,data,forward)
    type(DFTI_DESCRIPTOR), pointer, intent(in) :: plan
    complex(mytype), intent(inout) :: data(*)
    logical, intent(in) :: forward
    ! Sequence association flattens a contiguous Y slab or Z pencil without
    ! copying. Strides/batches above follow the installed 2DECOMP MKL backend.
    if (forward) then
       call check_fft(DftiComputeForward(plan,data))
    else
       call check_fft(DftiComputeBackward(plan,data))
    endif
  end subroutine transform
#endif

  subroutine scalar_fourier_finalize()
    if (.not.initialized) return
#ifdef X3D_SCALAR_MAPPING_MKL
    call check_fft(DftiFreeDescriptor(plan_y))
    call check_fft(DftiFreeDescriptor(plan_z))
#endif
    deallocate(sx,sy,sz,inverse_pivot,width)
    initialized=.false.
  end subroutine scalar_fourier_finalize

end module scalar_mapping_fourier
