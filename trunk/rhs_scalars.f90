subroutine rhs_scalars

  use m_openmpi
  use m_io
  use m_parameters
  use m_fields
  use m_work
  use x_fftw

  implicit none

  integer :: i, j, k, n
  real*8  :: rtmp1, rtmp2, wnum2, akmax

  ! converting them to the real space
  wrk(:,:,:,1:3) = fields(:,:,:,1:3)
  do n = 1,3
     call xFFT3d(-1,n)
  end do

  ! If we're not advancing scalars, return
  if (.not.int_scalars) return


  ! Do each scalar at a time

  ! Trying to keep the velocities in wrk1:3 intact because they
  ! are needed later 
  do n = 1, n_scalars

     wrk(:,:,:,0) = fields(:,:,:,3+n)
     call xFFT3d(-1,0)

     ! Products of the scalar and velocities
     do i = 1,3
        wrk(:,:,:,n+2+i) = - wrk(:,:,:,0) * wrk(:,:,:,i)
        call xFFT3d(1,n+2+i)
     end do

     ! kmax in real
     akmax = real(kmax,8)

     ! Assembling the RHS in wrk(:,:,:,3+n)
     do k = 1,nz
        do j = 1,ny
           do i = 1,nx+1,2

!  -------------- Edit 3
              ! If the dealiasing option is 2/3-rule (dealias=0) then we retain the modes
              ! inside the cube described by $| k_i | \leq  k_{max}$, $i=1,2,3$.
              ! The rest of the modes is purged


              ! Only bothering with the wavenumbers not greater than kmax
              wnum2 = akx(i)**2 + aky(k)**2 + akz(j)**2
              if (wnum2 .gt. real(kmax**2,8)) then

                 ! all the wavenumbers that are greater than kmax get zeroed out
                 wrk(i  ,j,k,3+n) = zip
                 wrk(i+1,j,k,3+n) = zip

              else
                 ! taking the convective term, multiply it by "i" 
                 ! (see how it's done in x_fftw.f90)
                 ! and adding the diffusion term

                 ! also using the fact that the waveunmbers for (i,j,k) are the same
                 ! as wavenumbers for (i+1,j,k)

                 ! i * (a + ib) + d = -b + ia + d
                 rtmp1 =   akx(i+1)*wrk(i+1,j,k,3+n) + aky(k)*wrk(i+1,j,k,4+n) + akz(j)*wrk(i+1,j,k,5+n)
                 rtmp2 =   akx(i  )*wrk(i  ,j,k,3+n) + aky(k)*wrk(i  ,j,k,4+n) + akz(j)*wrk(i  ,j,k,5+n)

                 ! also using the fact that the waveunmbers for (i,j,k) are the same
                 ! as wavenumbers for (i+1,j,k)

                 wrk(i  ,j,k,3+n) = - rtmp1 - pe(n) * wnum2*fields(i  ,j,k,3+n)
                 wrk(i+1,j,k,3+n) =   rtmp2 - pe(n) * wnum2*fields(i+1,j,k,3+n)

              end if
           end do
        end do
     end do

     ! Now adding the reaction part
     if (scalar_type(n).ge.100) then
        call add_reaction(n)
        call dealias_rhs(3+n)
     end if

  end do

  return
end subroutine rhs_scalars

!================================================================================
!================================================================================
subroutine add_reaction(n)

  use m_openmpi
  use m_io
  use m_parameters
  use m_fields
  use m_work
  use x_fftw

  implicit none

  integer :: n, rtype
  real*8  :: scmean, rrate

  ! reaction type
  rtype =  scalar_type(n)/100

  ! raction rate 
  rrate = reac_sc(n)

  select case (rtype)
  case (1)

     ! KPP reaction rate
     wrk(:,:,:,0) = rrate * (1.d0 - wrk(:,:,:,0)**2)

  case (2)

     ! symmetric bistable       
     wrk(:,:,:,0) = rrate * (1.d0 - wrk(:,:,:,0)**2) * wrk(:,:,:,0)

  case (3)

     ! self-adjusting bistable
     scmean = fields(1,1,1,3+n)/nxyz_all
     call MPI_BCAST(scmean, 1, MPI_REAL8, 0, MPI_COMM_TASK, mpi_err)

     wrk(:,:,:,0) = rrate * (1.d0 - wrk(:,:,:,0)**2) * &
          (wrk(:,:,:,0) - scmean)
  case default

     write(out,*) "Unknown reaction rate"
     call flush(out)

     stop

  end select


  ! FFT the reaction into the Fourier space
  call xFFT3d(1,0)

  ! Adding reaction to the RHS in wrk(:,:,:,3+n)

  wrk(:,:,:,3+n) = wrk(:,:,:,3+n) + wrk(:,:,:,0)

end subroutine add_reaction

!================================================================================

subroutine dealias_rhs(n)

  use m_io
  use m_parameters
  use m_work
  use x_fftw

  implicit none

  integer :: i, j, k, n
  real*8  :: wnum2, akmax

  akmax = real(kmax,8)

  do k = 1,nz
     do j = 1,ny
        do i = 1,nx+1,2

!!$           wnum2 = akx(i)**2 + aky(k)**2 + akz(j)**2
!!$           if (wnum2 .gt. real(kmax**2,8)) then

           if  (abs(akx(i)) .gt. akmax .or. &
                abs(aky(k)) .gt. akmax .or. &
                abs(akz(j)) .gt. akmax) then
              wrk(i  ,j,k,n) = zip
              wrk(i+1,j,k,n) = zip
           end if

        end do
     end do
  end do

  return

end subroutine dealias_rhs

!================================================================================
!================================================================================
!================================================================================
!================================================================================
subroutine test_rhs_scalars

  use m_openmpi
  use m_io
  use m_parameters
  use m_fields
  use m_work
  use x_fftw

  implicit none

  integer :: i,j,k, n
  real*8 :: a,b,c, x,y,z
  if (task.eq.'hydro') then

     a = 1.d0
     b = 5.d0
     c = 17.d0

     do k = 1,nz
        do j = 1,ny
           do i = 1,nx

              x = dx*real(i-1)
              y = dx*real(j-1)
              z = dx*real(myid*nz + k-1)

              wrk(i,j,k,1) = sin(a * x)
              wrk(i,j,k,2) = sin(b * y)
              wrk(i,j,k,3) = sin(c * z)
              wrk(i,j,k,4) = cos(a * x)
           end do
        end do
     end do


     do n = 1,4
        call xFFT3d(1,n)
        fields(:,:,:,n) = wrk(:,:,:,n)
     end do


     nu = .5d0

     call rhs_scalars

     print *,'got rhs'

     call xFFT3d(-1,4)


     do k = 1,nz
        do j = 1,ny
           do i = 1,nx

              x = dx*real(i-1)
              y = dx*real(j-1)
              z = dx*real(myid*nz + k-1)


              ! checking 
              wrk(i,j,k,0) = -a*cos(2.*a*x) - cos(a*x)*(b*cos(b*y) + c*cos(c*z) + nu*a**2)

           end do
        end do
     end do

!!$     tmp4(:,:,:) = wrk(1:nx,:,:,4)
!!$     fname = 'r1.arr'
!!$     call write_tmp4
!!$
!!$     tmp4(:,:,:) = wrk(1:nx,:,:,0)
!!$     fname = 'r0.arr'
!!$     call write_tmp4




     wrk(:,:,:,0) = abs(wrk(1:nx,:,:,0) - wrk(1:nx,:,:,4))

     print *,'Maximum error is ',maxval(wrk(1:nx,:,:,0))

!!$     tmp4(:,:,:) = wrk(1:nx,:,:,3) - wrk(1:nx,:,:,6)
!!$     fname = 'e3.arr'
!!$     call write_tmp4
!!$


  end if
  return
end subroutine test_rhs_scalars
