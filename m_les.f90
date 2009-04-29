!================================================================================
!  m_les - the module that contains all subroutines and variables that are
!  needed to introduce Large Eddy Simulation (LES) into the code.
!
!  The behaviour of the module is governed by the variable "les_mode" from the
!  module m_parameters.f90
!
!  Time-stamp: <2009-04-22 16:34:24 (chumakov)>
!================================================================================
module m_les

  use m_openmpi
  use m_io
  use m_parameters
  use m_fields
  use m_work
  implicit none


!================================================================================
!  Arrays and constants
!================================================================================

  ! indicator, whether we use LES at all
  logical :: les = .false.

  ! model switch "les_mode" is contained by m_parameters
  ! integer(kind=4) :: les_model

  ! array for turbulent viscosity
  real(kind=8), allocatable :: turb_visc(:,:,:)

  ! Smagorinsky constant
  real(kind=8) :: c_smag = 0.18

  ! test filter width
  real(kind=8) :: les_delta

  ! array with the test filter in it
  real(kind=8), allocatable :: filter_g(:,:,:)

  ! model indicator for the output
  character*3 :: les_model_name = '   '

!================================================================================
!                            SUBROUTINES
!================================================================================
contains
!================================================================================
!  Allocation of LES arrays
!================================================================================
  subroutine m_les_init
    implicit none

    real*8, allocatable :: sctmp(:)

    ! if les_model=0, do not initialize anything and return
    if (les_model==0) return

    write(out,*) 'Initializing LES...'
    call flush(out)

    ! depending on the value of les_model, initialize different things
    les = .true.
    ! initialize the filter width to be equal to the grid spaxing
    les_delta = dx
    ! initializeing stuff based on the model switch "les_model"
    select case (les_model)
    case(1)
       ! Smagorinsky model
       write(out,*) ' - Smagorinsky model: initializing the eddy viscosity'
       call flush(out)
       allocate(turb_visc(nx,ny,nz),stat=ierr)
       if (ierr/=0) then
          write(out,*) 'Cannot allocate the turbulent viscosity array'
          call flush(out)
          call my_exit(-1)
       end if
       turb_visc = 0.0d0

       n_les = 0
       les_model_name = " SM"

    case(2)
       ! Dynamic Localization model
       write(out,*) ' - DL model: initializing the eddy viscosity and adding extra transport equation'
       call flush(out)
       allocate(turb_visc(nx,ny,nz),stat=ierr)
       if (ierr/=0) then
          write(out,*) 'Cannot allocate the turbulent viscosity array'
          call flush(out)
          call my_exit(-1)
       end if
       turb_visc = 0.0d0

       n_les = 1
       les_model_name = "DLM"

    case default
       write(out,*) 'M_LES_INIT: invalid value of les_model:',les_model
       call flush(out)
       call my_exit(-1)
    end select

    write(out,*) "n_les = ", n_les

    additional_scalars: if (n_les .gt. 0) then
       write(out,*) "Adding elements to arrays PE and SC for the LES-related scalars..."
       call flush(out)
       allocate(sctmp(1:n_scalars+n_les), stat=ierr)
       sctmp(1:n_scalars) = sc(1:n_scalars)
       sctmp(n_scalars+1:n_scalars+n_les) = one
       if (allocated(sc)) deallocate(sc)
       allocate(sc(1:n_scalars+n_les))
       sc = sctmp
       if (allocated(pe)) deallocate(pe)
       allocate(pe(1:n_scalars+n_les))
       pe = nu / sc
       deallocate(sctmp)
       write(out,*) " ...done."
       call flush(out)
    end if additional_scalars


    write(out,*) 'initialized LES.'
    call flush(out)

    return
  end subroutine m_les_init


!================================================================================
!  initialization of the LES arrays - part 2
!  definition of the arrays
!
!  called after the restart
!================================================================================
  subroutine m_les_begin

    implicit none

    ! if les_model=0, do not initialize anything and return
    ! also don't do anything if the model is Smagorinsky model - it does not
    ! need any additional initialization beside the array allocation which has
    ! been already done.
    if (les_model<=1) return

    write(out,*) "M_LES_BEGIN..."
    call flush(out)

    select case(les_model)

    case(2)
       write(out,*) "-- DLM model"
       write(out,*) "-- Initializing k_sgs"
       call flush(out)

       call m_les_dlm_k_init


    case default
       write(out,*) "M_LES_BEGIN: invalid value of les_model: ",les_model
       call flush(out)
       call my_exit(-1)
    end select


    return
  end subroutine m_les_begin

!================================================================================
!  Adding LES sources to the RHS of velocities
!================================================================================
  subroutine les_rhs_velocity

    implicit none

!!$    write(out,*) "Subroutine LES_RHS_VELOCITY"
!!$    call flush(out)

    select case (les_model)
    case(1:2)
!!$       write(out,*) ' Turbulent viscosity (Smag or DLM)'
!!$       call flush(out)
       call les_rhsv_turb_visc
    case default
       write(out,*) 'LES_RHS_VELOCITY: invalid value of les_model:',les_model
       call flush(out)
       call my_exit(-1)
    end select

    return
  end subroutine les_rhs_velocity

!================================================================================
!================================================================================
!  Adding LES sources to the RHS of scalars
!================================================================================
  subroutine les_rhs_scalars
    implicit none

    select case (les_model)
    case(1)
       call les_rhss_turb_visc
    case(2)
       call m_les_k_src_dlm
       call m_les_k_diss_algebraic
       call les_rhss_turb_visc
    case default
       write(out,*) 'LES_RHS_SCALARS: invalid value of les_model:',les_model
       call flush(out)
       call my_exit(-1)
    end select
    return
  end subroutine les_rhs_scalars

!================================================================================
!    calculating velocity sources and adding them to the RHS's (wrk1...3)
!
!    case when SGS stress tau_ij is modeled using turbulent viscosity
!    the turb. viscosity is supposed to be in the array turb_visc(nx,ny,nz) 
!================================================================================
  subroutine les_rhsv_turb_visc

    use x_fftw
    implicit none


    integer :: i, j, k, n
    real(kind=8) :: rtmp

    ! due to memory constraints we have only three work arrays wrk4..6,
    ! because the first three wrk1..3 contain already comptued velocity RHS's.

    ! Calculating S_11, S_12, S_13
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx + 1, 2
             ! S_11, du/dx
             wrk(i  ,j,k,4) = - akx(i+1) * fields(i+1,j,k,1)
             wrk(i+1,j,k,4) =   akx(i  ) * fields(i  ,j,k,1)

             ! S_12, 0.5 (du/dy + dv/dx)
             wrk(i  ,j,k,5) = -half * ( aky(k) * fields(i+1,j,k,1) + akx(i+1) * fields(i+1,j,k,2) )
             wrk(i+1,j,k,5) =  half * ( aky(k) * fields(i  ,j,k,1) + akx(i  ) * fields(i  ,j,k,2) )

             ! S_13, 0.5 (du/dz + dw/dx)
             wrk(i  ,j,k,6) = -half * ( akz(j) * fields(i+1,j,k,1) + akx(i+1) * fields(i+1,j,k,3) )
             wrk(i+1,j,k,6) =  half * ( akz(j) * fields(i  ,j,k,1) + akx(i  ) * fields(i  ,j,k,3) )
          end do
       end do
    end do

    ! Multiplying them by  -2 * turbulent viscosity to get tau_11, tau_12, tau_13
    do n = 4,6
       call xFFT3d(-1,n)
       wrk(1:nx,1:ny,1:nz,n) = - two * turb_visc(1:nx,1:ny,1:nz) * wrk(1:nx,1:ny,1:nz,n) 
       call xFFT3d(1,n)
    end do

    ! Taking d/dx tau_11,  d/dy tau_12, d/dz tau_13 and subtracting from the current RHS
    ! note the sign reversal (-/+) because we subtract this from RHS
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx+1, 2

             ! Cutting off any wave modes that can introduce aliasing into the velocities
             ! This "dealiasing" is done in the most restrictive manner: only Fourier modes
             ! that have ialias=0 are added
             if (ialias(i,j,k) .eq. 0) then

                wrk(i  ,j,k,1) = wrk(i  ,j,k,1) + akx(i+1)*wrk(i+1,j,k,4) + aky(k)*wrk(i+1,j,k,5) + akz(j)*wrk(i+1,j,k,6)
                wrk(i+1,j,k,1) = wrk(i+1,j,k,1) - akx(i  )*wrk(i  ,j,k,4) - aky(k)*wrk(i  ,j,k,5) - akz(j)*wrk(i  ,j,k,6)

                wrk(i  ,j,k,2) = wrk(i  ,j,k,2) + akx(i+1)*wrk(i+1,j,k,5)
                wrk(i+1,j,k,2) = wrk(i+1,j,k,2) - akx(i  )*wrk(i  ,j,k,5)

                wrk(i  ,j,k,3) = wrk(i  ,j,k,3) + akx(i+1)*wrk(i+1,j,k,6)
                wrk(i+1,j,k,3) = wrk(i+1,j,k,3) - akx(i  )*wrk(i  ,j,k,6)

             end if

          end do
       end do
    end do


    ! Calculating S_22, S_23, S_33
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx+1, 2

             ! S_22, dv/dy
             wrk(i  ,j,k,4) = - aky(k) * fields(i+1,j,k,2)
             wrk(i+1,j,k,4) =   aky(k) * fields(i  ,j,k,2)

             ! S_23, 0.5 (dv/dz + dw/dy)
             wrk(i  ,j,k,5) = - half * ( akz(j) * fields(i+1,j,k,2) + aky(k) * fields(i+1,j,k,3) )
             wrk(i+1,j,k,5) =   half * ( akz(j) * fields(i  ,j,k,2) + aky(k) * fields(i  ,j,k,3) )

             ! S_33, de/dz
             wrk(i  ,j,k,6) = - akz(j) * fields(i+1,j,k,3)
             wrk(i+1,j,k,6) =   akz(j) * fields(i  ,j,k,3)
          end do
       end do
    end do

    ! Multiplying them by -2 * turbulent viscosity to get tau_22, tau_23, tau_33
    do n = 4,6
       call xFFT3d(-1,n)
       wrk(1:nx,1:ny,1:nz,n) = - two * turb_visc(1:nx,1:ny,1:nz) * wrk(1:nx,1:ny,1:nz,n) 
       call xFFT3d(1,n)
    end do

    ! Taking
    ! d/dy tau_22, d/dz tau_23
    ! d/dy tau_23, d/dz tau_33 
    ! and subtracting from the current RHS for v and w
    ! note the sign reversal (-/+) because we subtract this from RHS
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx+1, 2

             ! Cutting off any wave modes that can introduce aliasing into the velocities
             ! This "dealiasing" is done in the most restrictive manner: only Fourier modes
             ! that have ialias=0 are added
             if (ialias(i,j,k) .eq. 0) then
                wrk(i  ,j,k,2) = wrk(i  ,j,k,2) + aky(k)*wrk(i+1,j,k,4) + akz(j)*wrk(i+1,j,k,5)
                wrk(i+1,j,k,2) = wrk(i+1,j,k,2) - aky(k)*wrk(i  ,j,k,4) - akz(j)*wrk(i  ,j,k,5)

                wrk(i  ,j,k,3) = wrk(i  ,j,k,3) + aky(k)*wrk(i+1,j,k,5) + akz(j)*wrk(i+1,j,k,6)
                wrk(i+1,j,k,3) = wrk(i+1,j,k,3) - aky(k)*wrk(i  ,j,k,5) - akz(j)*wrk(i  ,j,k,6)
             end if

          end do
       end do
    end do

    return
  end subroutine les_rhsv_turb_visc

!================================================================================
!    calculating LES sources for scalars and adding them to the RHS's 
!    wrk4...3+n_scalars+n_les
!
!    case when SGS stress tau_ij is modeled using turbulent viscosity
!    the turb. viscosity is supposed to be in the array turb_visc(nx,ny,nz) 
!================================================================================
  subroutine les_rhss_turb_visc

    use x_fftw
    implicit none

    integer :: i, j, k, n, tmp1, tmp2
    character :: dir
    real(kind=8) :: rtmp

    ! have two arrays available as work arrays: wrk 3+n_scalars+n_les+1 and +2
    tmp1 = 3 + n_scalars + n_les + 1
    tmp2 = 3 + n_scalars + n_les + 2

    ! For every scalar (passive or LES active scuch as k_sgs), add turbulent
    ! viscosity to the RHS
    les_rhs_all_scalars: do n = 4, 3 + n_scalars + n_les

!!$       write(out,*) "les_rhss_turb_visc: doing field #",n
!!$       call flush(out)

       ! computing the second derivative, multiplying it by turb_visc
       ! and adding to the RHS (that is contained in wrk(n))

       wrk(:,:,:,tmp1) = fields(:,:,:,n)
       wrk(:,:,:,0) = zip

       directions: do i = 1,3
          if (i.eq.1) dir = 'x'
          if (i.eq.2) dir = 'y'
          if (i.eq.3) dir = 'z'

          ! taking the first derivatite, multiplying by the turb_visc
          ! then taking another derivative and adding the result to wrk0
          call x_derivative(tmp1, dir, tmp2)
          call xFFT3d(-1,tmp2)
          wrk(1:nx, 1:ny, 1:nz, tmp2) = wrk(1:nx, 1:ny, 1:nz, tmp2) * turb_visc(1:nx, 1:ny, 1:nz)
          call xFFT3d(1,tmp2)
          call x_derivative(tmp2, dir, tmp2)

          ! following Yoshizawa and Horiuti (1985), the viscosity in the scalar equation
          ! is twice the viscosity used in the production of k_sgs
          ! (Journal of Phys. Soc. Japan, V.54 N.8, pp.2834-2839)
          wrk(1:nx, 1:ny, 1:nz, tmp2) = two * wrk(1:nx, 1:ny, 1:nz, tmp2)

          wrk(:,:,:,0) = wrk(:,:,:,0) + wrk(:,:,:,tmp2)
       end do directions

       ! adding the d/dx ( nu_t d phi/dx) to the RHS for the scalar (field #n)
       ! only adding the Fourier modes that are not producing any aliasing
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx+2
                if (ialias(i,j,k) .eq. 0) wrk(i,j,k,n) = wrk(i,j,k,n) + wrk(i,j,k,0)
             end do
          end do
       end do

    end do les_rhs_all_scalars

    return
  end subroutine les_rhss_turb_visc


!================================================================================
!================================================================================
!  Calculation of turbulent viscosity turb_visc(:,:,:)
!================================================================================
  subroutine les_get_turb_visc

    implicit none

    select case (les_model)
    case(1)
       call les_get_turb_visc_smag

    case(2)
       call les_get_turb_visc_dlm

    case default
       write(out,*) "LES_GET_TURB_VISC: les_model: ", les_model
       write(out,*) "LES_GET_TURB_VISC: Not calculating turb_visc"
       call flush(out)
       call my_exit(-1)
    end select

    return
  end subroutine les_get_turb_visc


!================================================================================
!  Calculation of turbulent viscosity turb_visc(:,:,:) - Smagorinsky model
!================================================================================
!================================================================================
  subroutine les_get_turb_visc_smag

    use x_fftw
    implicit none

    integer :: i, j, k, n
    real*8 :: c_smag = 0.18_8

    ! due to memory constraints we have only three work arrays wrk4..6,
    ! because the first three wrk1..3 contain already comptued velocity RHS's.

    ! Calculating S_11, S_12, S_13
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx + 1, 2
             ! S_11, du/dx
             wrk(i  ,j,k,4) = - akx(i+1) * fields(i+1,j,k,1)
             wrk(i+1,j,k,4) =   akx(i  ) * fields(i  ,j,k,1)

             ! S_12, 0.5 (du/dy + dv/dx)
             wrk(i  ,j,k,5) = -half * ( aky(k) * fields(i+1,j,k,1) + akx(i+1) * fields(i+1,j,k,2) )
             wrk(i+1,j,k,5) =  half * ( aky(k) * fields(i  ,j,k,1) + akx(i  ) * fields(i  ,j,k,2) )

             ! S_13, 0.5 (du/dz + dw/dx)
             wrk(i  ,j,k,6) = -half * ( akz(j) * fields(i+1,j,k,1) + akx(i+1) * fields(i+1,j,k,3) )
             wrk(i+1,j,k,6) =  half * ( akz(j) * fields(i  ,j,k,1) + akx(i  ) * fields(i  ,j,k,3) )
          end do
       end do
    end do

    ! Converting them to real space and adding to turb_visc(:,:,:)
    do n = 4,6
       call xFFT3d(-1,n)
    end do
    turb_visc(1:nx,1:ny,1:nz) =                                   wrk(1:nx,1:ny,1:nz,4)**2
    turb_visc(1:nx,1:ny,1:nz) = turb_visc(1:nx,1:ny,1:nz) + two * wrk(1:nx,1:ny,1:nz,5)**2
    turb_visc(1:nx,1:ny,1:nz) = turb_visc(1:nx,1:ny,1:nz) + two * wrk(1:nx,1:ny,1:nz,6)**2

    ! Calculating S_22, S_23, S_33
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx+1, 2

             ! S_22, dv/dy
             wrk(i  ,j,k,4) = - aky(k) * fields(i+1,j,k,2)
             wrk(i+1,j,k,4) =   aky(k) * fields(i  ,j,k,2)

             ! S_23, 0.5 (dv/dz + dw/dy)
             wrk(i  ,j,k,5) = - half * ( akz(j) * fields(i+1,j,k,2) + aky(k) * fields(i+1,j,k,3) )
             wrk(i+1,j,k,5) =   half * ( akz(j) * fields(i  ,j,k,2) + aky(k) * fields(i  ,j,k,3) )

             ! S_33, dw/dz
             wrk(i  ,j,k,6) = - akz(j) * fields(i+1,j,k,3)
             wrk(i+1,j,k,6) =   akz(j) * fields(i  ,j,k,3)
          end do
       end do
    end do

    ! Converting them to real space and adding to turb_visc(:,:,:)
    do n = 4,6
       call xFFT3d(-1,n)
    end do
    turb_visc(1:nx,1:ny,1:nz) = turb_visc(1:nx,1:ny,1:nz) +       wrk(1:nx,1:ny,1:nz,4)**2
    turb_visc(1:nx,1:ny,1:nz) = turb_visc(1:nx,1:ny,1:nz) + two * wrk(1:nx,1:ny,1:nz,5)**2
    turb_visc(1:nx,1:ny,1:nz) = turb_visc(1:nx,1:ny,1:nz) +       wrk(1:nx,1:ny,1:nz,6)**2
    ! now turb_visc contains S_ij S_ij

    ! Finishing up the turbulent viscosiy
    ! making it (C_s Delta)^2 |S|, where |S| = sqrt(2 S_{ij} S_{ij})
    turb_visc = sqrt( two * turb_visc)
    turb_visc  = turb_visc * (c_smag * les_delta)**2

    return
  end subroutine les_get_turb_visc_smag

!================================================================================
!================================================================================
!  Calculation of turbulent viscosity turb_visc(:,:,:) - DLM model
!================================================================================
  subroutine les_get_turb_visc_dlm

    use x_fftw
    implicit none

!!$    real*8 :: C_k = 0.05d0 ! This is take from Yoshizawa and Horiuti (1985)
    real*8 :: C_k = 0.1d0 ! This is what works for this code.  Dunno why...
    real*8 :: sctmp, sctmp1

!!$    write(out,*) "Calculating turbulent viscosity using DLM model"
!!$    call flush(out)

    wrk(:,:,:,0) = fields(:,:,:,3+n_scalars+1)
    call xFFT3d(-1,0)

    ! C alculating the minimum value of k_sgs.  If it's less that zero, clip it
    ! at zero.
    sctmp1 = minval(wrk(1:nx,1:ny,1:nz,0))
    count = 1
    call MPI_REDUCE(sctmp1,sctmp,count,MPI_REAL8,MPI_MIN,0,MPI_COMM_TASK,mpi_err)
    call MPI_BCAST(sctmp,count,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)
    if (sctmp.lt.zip) then
       ! write(out,*) 'LES_GET_TURB_VISC_DLM: minval of k is less than 0:',sctmp
       ! call flush(out)
       ! call my_exit(-1)
       wrk(:,:,:,0) = max(wrk(:,:,:,0), zip)       
    end if

    turb_visc(1:nx,1:ny,1:nz) = C_k * les_delta * sqrt(wrk(1:nx,1:ny,1:nz,0))

    ! If the minimum of k_sgs is less than zero, we want to clip it
    ! mercilessly (already done in wrk0) then put it back in Fourier space
    if (sctmp.lt.zip) then
       call xFFT3d(1,0)
       fields(:,:,:,3+n_scalars+1) = wrk(:,:,:,0)
    end if

!!$    write(out,*) "      turbulent viscosity calculated"
!!$    call flush(out)

    return
  end subroutine les_get_turb_visc_dlm


!================================================================================
!  Subroutine that initializes k_sgs for the Dynamic Localization Model (DLM)
!================================================================================
  subroutine m_les_dlm_k_init

    use x_fftw
    implicit none

    integer :: i, j, k, n_k

    write(out,*) "m_les_dlm_k_init: initializing k_sgs"
    call flush(out)

    n_k = 3 + n_scalars + 1
    fields(:,:,:,n_k) = zip

    ! initializing it as a constant
    write(out,*) "m_les_dlm_k_init: initialized k=0.2"
    call flush(out)
    fields(:,:,:,n_k) = 0.2d0
    call xFFT3d_fields(1,n_k)
    return


    do i = 1, 3
       wrk(:,:,:,1) = fields(:,:,:,i)
       call x_derivative(1, 'x', 2)
       call x_derivative(1, 'y', 3)
       call x_derivative(1, 'z', 4)

       call xFFT3d(-1,2)
       call xFFT3d(-1,3)
       call xFFT3d(-1,4)

       fields(:,:,:,n_k) = fields(:,:,:,n_k) + wrk(:,:,:,2)**2
       fields(:,:,:,n_k) = fields(:,:,:,n_k) + wrk(:,:,:,3)**2
       fields(:,:,:,n_k) = fields(:,:,:,n_k) + wrk(:,:,:,4)**2
    end do

    ! put k_sgs in real space into fields(:,:,:,n_k)
    fields(:,:,:,n_k) = fields(:,:,:,n_k) * les_delta**2 / 12.d0

    write(out,*) "m_les_dlm_k_init: minval(k)=",minval(fields(1:nx,1:ny,1:nz,n_k))
    call flush(out)


    ! convert to Fourier space and dealias
    call xFFT3d_fields(1,n_k)

    do k = 1, nz
       do j = 1, ny
          do i = 1, nx+2
             if (ialias(i,j,k) .gt. 0) fields(i,j,k,n_k) = zip
          end do
       end do
    end do

    wrk(:,:,:,0) = fields(:,:,:,n_k)
    call xFFT3d(-1,0)
    write(out,*) "m_les_dlm_k_init: minval(k)=",minval(wrk(1:nx,1:ny,1:nz,0))
    call flush(out)


    return
  end subroutine m_les_dlm_k_init

!================================================================================
!================================================================================
!  Subroutine that calculates the source for k_sgs:   - tau_{ij} S_{ij} 
!  for the case of Dynamic Localization Model (DLM)
!================================================================================
  subroutine m_les_k_src_dlm

    use x_fftw
    implicit none
    integer :: n1, n2, k_n, i, j, k

    ! The source itself is nothing but 2 * nu_t * S_{ij} * S_{ij}
    ! So the main hassle is to calculate S_{ij}, since nu_t is available from
    ! the array turb_visc.

    ! have two arrays available as work arrays: wrk 3+n_scalars+n_les+1 and +2
    ! also have wrk0 in which we will assemble the source at the end
    n1 = 3 + n_scalars + n_les + 1
    n2 = 3 + n_scalars + n_les + 2

    wrk(:,:,:,0) = zip
    wrk(:,:,:,n1) = zip
    wrk(:,:,:,n2) = zip

    ! calculating the S_{ij} S_{ij}.  Note that when calculating derivatives,
    ! we only process those Fourier models that won't introduce aliasing when
    ! the wuantity is squared.  These modes are given by ialias(i,j,k)=0

    ! calculating S_11, S_12
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx + 1, 2
             if (ialias(i,j,k).eq.0) then
                ! S_11, du/dx
                wrk(i  ,j,k,n1) = - akx(i+1) * fields(i+1,j,k,1)
                wrk(i+1,j,k,n1) =   akx(i  ) * fields(i  ,j,k,1)
                ! S_12, 0.5 (du/dy + dv/dx)
                wrk(i  ,j,k,n2) = -half * ( aky(k) * fields(i+1,j,k,1) + akx(i+1) * fields(i+1,j,k,2) )
                wrk(i+1,j,k,n2) =  half * ( aky(k) * fields(i  ,j,k,1) + akx(i  ) * fields(i  ,j,k,2) )
             end if
          end do
       end do
    end do
    ! converting to real space, squaring and adding to wrk0
    call xFFT3d(-1,n1);  call xFFT3d(-1,n2); 
    wrk(:,:,:,0) = wrk(:,:,:,0) + wrk(:,:,:,n1)**2 + two*wrk(:,:,:,n2)**2

    ! calculating S_13, S_22
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx + 1, 2
             if(ialias(i,j,k).eq.0) then
                ! S_13, 0.5 (du/dz + dw/dx)
                wrk(i  ,j,k,n1) = -half * ( akz(j) * fields(i+1,j,k,1) + akx(i+1) * fields(i+1,j,k,3) )
                wrk(i+1,j,k,n1) =  half * ( akz(j) * fields(i  ,j,k,1) + akx(i  ) * fields(i  ,j,k,3) )
                ! S_22, dv/dy
                wrk(i  ,j,k,n2) = - aky(k) * fields(i+1,j,k,2)
                wrk(i+1,j,k,n2) =   aky(k) * fields(i  ,j,k,2)
             end if
          end do
       end do
    end do
    ! converting to real space, squaring and adding to wrk0
    call xFFT3d(-1,n1);  call xFFT3d(-1,n2); 
    wrk(:,:,:,0) = wrk(:,:,:,0) + two*wrk(:,:,:,n1)**2 + wrk(:,:,:,n2)**2

    ! calculating S_23, S_33
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx + 1, 2
             if(ialias(i,j,k).eq.0) then
                ! S_23, 0.5 (dv/dz + dw/dy)
                wrk(i  ,j,k,n1) = - half * ( akz(j) * fields(i+1,j,k,2) + aky(k) * fields(i+1,j,k,3) )
                wrk(i+1,j,k,n1) =   half * ( akz(j) * fields(i  ,j,k,2) + aky(k) * fields(i  ,j,k,3) )
                ! S_33, de/dz
                wrk(i  ,j,k,n2) = - akz(j) * fields(i+1,j,k,3)
                wrk(i+1,j,k,n2) =   akz(j) * fields(i  ,j,k,3)
             end if
          end do
       end do
    end do
    ! converting to real space, squaring and adding to wrk0
    call xFFT3d(-1,n1);  call xFFT3d(-1,n2); 
    wrk(:,:,:,0) = wrk(:,:,:,0) + two*wrk(:,:,:,n1)**2 + wrk(:,:,:,n2)**2

    ! at this point wrk0 contains S_{ij} S_{ij} in real space.
    ! need to multiply by two and multiply by turb_visc to get the source
    ! (the energy transfer term)
    wrk(:,:,:,0) = two * wrk(:,:,:,0)
    wrk(1:nx,1:ny,1:nz,0) = wrk(1:nx,1:ny,1:nz,0) * turb_visc(1:nx,1:ny,1:nz)

    ! convert the transfer term to Fourier space
    call xFFT3d(1,0)

    ! adding this energy transfer term to the RHS for k_sgs
    ! the RHS for k_sgs is supposed to be in wrk(3+n_scalars+1)
    k_n = 3 + n_scalars + 1

    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             if (ialias(i,j,k).eq.0) wrk(i,j,k,k_n) = wrk(i,j,k,k_n) + wrk(i,j,k,0)
          end do
       end do
    end do
!!$    wrk(:,:,:,k_n) = wrk(:,:,:,k_n) + wrk(:,:,:,0)

  end subroutine m_les_k_src_dlm

!================================================================================
!================================================================================
!  Dissipation term in k-equation: simple algebraic model k^(3/2)/Delta
!================================================================================
  subroutine m_les_k_diss_algebraic

    use x_fftw
    implicit none
    integer :: n_k, i, j, k
    real*8 :: sctmp, sctmp1

    ! the "field number" for k_sgs
    n_k = 3 + n_scalars + 1

    ! get the SGS kinetic energy in x-space
    wrk(:,:,:,0) = fields(:,:,:,n_k)
    
   ! first zero the modes that can produce aliasing 
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx+2
             if (ialias(i,j,k).gt.0) wrk(i,j,k,0) = zip
          end do
       end do
    end do
    ! then convert to x-space
    call xFFT3d(-1,0)

    ! check the minimum value of k
    sctmp1 = minval(wrk(1:nx,1:ny,1:nz,0))
    count = 1
    call MPI_REDUCE(sctmp1,sctmp,count,MPI_REAL8,MPI_MIN,0,MPI_COMM_TASK,mpi_err)
    call MPI_BCAST(sctmp,count,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)
!!$    if (sctmp.lt.zip) then
!!$       write(out,*) itime,"Minimum value of k is ", sctmp
!!$       call flush(out)
!!$    end if

    ! calculating the dissipation rate
    wrk(:,:,:,0) = max(zip, wrk(:,:,:,0))
    wrk(:,:,:,0) = wrk(:,:,:,0)**1.5d0 / les_delta
    call xFFT3d(1,0)

    ! subtracting the dissipation rate from the RHS for k
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             if (ialias(i,j,k).eq.0) wrk(i,j,k,n_k) = wrk(i,j,k,n_k) - wrk(i,j,k,0)
          end do
       end do
    end do
!!$    wrk(:,:,:,n_k) = wrk(:,:,:,n_k) - wrk(:,:,:,0)

    return
  end subroutine m_les_k_diss_algebraic

!================================================================================
!================================================================================
end module m_les
