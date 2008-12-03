subroutine rhs_velocity

  use m_openmpi
  use m_io
  use m_parameters
  use m_fields
  use m_work
  use x_fftw

  implicit none

  integer :: i, j, k, n
  real*8  :: t1(6), rtmp, wnum2




  ! The IFFT of velocities has been done earlier in rhs_scalars
  ! the velocities were kept in wrk1...wrk3, intact.

!!$  ! putting the velocity field in the wrk array
!!$  wrk(:,:,:,1:3) = fields(:,:,:,1:3)
!!$  ! performing IFFT to convert them to the X-space
!!$  call xFFT3d(-1,1)
!!$  call xFFT3d(-1,2)
!!$  call xFFT3d(-1,3)


!-------------------------------------------------------------------------
  ! getting the Courant number (on the master process only)
  wrk(:,:,:,4) = abs(wrk(:,:,:,1)) + abs(wrk(:,:,:,2)) + abs(wrk(:,:,:,3))
  rtmp = maxval(wrk(1:nx,:,:,4))
  call MPI_REDUCE(rtmp,courant,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_TASK,mpi_err)
  if (variable_dt) then
     count = 1
     call MPI_BCAST(courant,count,MPI_REAL8,0,MPI_COMM_TASK,mpi_err)
  end if
  courant = courant * dt / dx

!-------------------------------------------------------------------------



  ! getting all 6 products of velocities
  do k = 1,nz
     do j = 1,ny
        do i = 1,nx
           t1(1) = wrk(i,j,k,1) * wrk(i,j,k,1)
           t1(2) = wrk(i,j,k,1) * wrk(i,j,k,2)
           t1(3) = wrk(i,j,k,1) * wrk(i,j,k,3)
           t1(4) = wrk(i,j,k,2) * wrk(i,j,k,2)
           t1(5) = wrk(i,j,k,2) * wrk(i,j,k,3)
           t1(6) = wrk(i,j,k,3) * wrk(i,j,k,3)

           do n = 1,6
              wrk(i,j,k,n) = t1(n)
           end do

        end do
     end do
  end do


  ! converting the products to the Fourier space
  do n = 1,6
     call xFFT3d(1,n)
  end do

  ! Building the RHS.  
  ! First, put into wrk arrays the convectove terms (that will be multiplyed by "i"
  ! later) and the factor that corresponds to the diffusion

  ! Do not forget that in Fourier space the indicies are (ix, iz, iy)
  do k = 1,nz
     do j = 1,ny
        do i = 1,nx+2

           t1(1) = - ( akx(i) * wrk(i,j,k,1) + aky(k) * wrk(i,j,k,2) + akz(j) * wrk(i,j,k,3) )
           t1(2) = - ( akx(i) * wrk(i,j,k,2) + aky(k) * wrk(i,j,k,4) + akz(j) * wrk(i,j,k,5) )
           t1(3) = - ( akx(i) * wrk(i,j,k,3) + aky(k) * wrk(i,j,k,5) + akz(j) * wrk(i,j,k,6) )


           t1(4) = - nu * ( akx(i)**2 + aky(k)**2 + akz(j)**2 ) 

           do n = 1,4
              wrk(i,j,k,n) = t1(n)
           end do

        end do
     end do
  end do

  ! now take the actual fields from fields(:,:,:,:) and calculate the RHSs

  ! at this moment the contains of wrk(:,:,:,1:3) are the convective terms in the RHS
  ! which are not yet multiplied by "i"
  ! wrk(:,:,:,4) contains the Laplace operator in Fourier space.  To get the diffusion term
  ! we need to take wrk(:,:,:,4) and multiply it by the velocity

  t1(5) = -nu * real(kmax,8)**2

  do k = 1,nz
     do j = 1,ny
        do i = 1,nx+1,2

           ! Only bothering with wavenumbers less than or equal to kmax
!!$           wnum2 = akx(i)**2 + aky(k)**2 + akz(j)**2
!!$           if (wnum2 .le. real(kmax**2,8)) then

           ! saving a little time by comparing two pre-computed numbers
           if ( t1(5) .le. wrk(i,j,k,4) ) then

              ! RHS for u, v and w
              do n = 1,3
                 ! taking the convective term, multiply it by "i" 
                 ! (see how it's done in x_fftw.f90)
                 ! and adding the diffusion term
                 rtmp =           - wrk(i+1,j,k,n) + wrk(i  ,j,k,4) * fields(i  ,j,k,n)
                 wrk(i+1,j,k,n) =   wrk(i  ,j,k,n) + wrk(i+1,j,k,4) * fields(i+1,j,k,n)
                 wrk(i  ,j,k,n) = rtmp
              end do

           else
              ! dealiasing: setting the Fourier components with magnitude higher than kmax to zero
              wrk(i  ,j,k,1:3) = zip
              wrk(i+1,j,k,1:3) = zip
           end if

        end do
     end do
  end do


  return
end subroutine rhs_velocity


!================================================================================
!================================================================================
!================================================================================
!================================================================================
!================================================================================
!================================================================================
!================================================================================

subroutine test_rhs_velocity  

  use m_openmpi
  use m_io
  use m_parameters
  use m_fields
  use m_work
  use x_fftw

  implicit none

  integer :: i,j,k, n
  real*8 :: a,b,c, x,y,z


  ! defining very particular velocities so the RHS can be computed analytically

  if (task.eq.'hydro') then


     write(out,*) 'inside.'
     call flush(out)


     a = 1.d0
     b = 1.d0
     c = 1.d0

     do k = 1,nz
        do j = 1,ny
           do i = 1,nx

              x = dx*real(i-1)
              y = dx*real(j-1)
              z = dx*real(myid*nz + k-1)

              wrk(i,j,k,1) = sin(a * x)
              wrk(i,j,k,2) = sin(b * y)
              wrk(i,j,k,3) = sin(c * z)
           end do
        end do
     end do

     write(out,*) 'did work'
     call flush(out)



     do n = 1,3
        call xFFT3d(1,n)
        fields(:,:,:,n) = wrk(:,:,:,n)
     end do


     write(out,*) 'did FFTs'
     call flush(out)

     nu = 0.d0

     call rhs_velocity

     write(out,*) 'got rhs'
     call flush(out)

     do n = 1,3
        call xFFT3d(-1,n)
     end do

     write(out,*) 'did FFTs'
     call flush(out)



     do k = 1,nz
        do j = 1,ny
           do i = 1,nx

              x = dx*real(i-1)
              y = dx*real(j-1)
              z = dx*real(myid*nz + k-1)

              ! checking u
              wrk(i,j,k,4) = -sin(a*x) * ( two*a*cos(a*x) + b*cos(b*y) + c*cos(c*z) + nu*a**2)

              ! checking v
              wrk(i,j,k,5) = -sin(b*y) * ( two*b*cos(b*y) + a*cos(a*x) + c*cos(c*z) + nu*b**2)

              ! checking w
              wrk(i,j,k,6) = -sin(c*z) * ( two*c*cos(c*z) + b*cos(b*y) + a*cos(a*x) + nu*c**2)


           end do
        end do
     end do

!!$    do k = 1,nz
!!$      write(out,"(3e15.6)") wrk(1,1,k,3),wrk(1,1,k,5),wrk(1,1,k,4)
!!$    end do

     wrk(:,:,:,0) = &
          abs(wrk(:,:,:,1) - wrk(:,:,:,4)) + &
          abs(wrk(:,:,:,2) - wrk(:,:,:,5)) + &
          abs(wrk(:,:,:,3) - wrk(:,:,:,6))

     print *,'Maximum error is ',maxval(wrk(1:nx,:,:,0))


     tmp4(:,:,:) = wrk(1:nx,:,:,1) - wrk(1:nx,:,:,4)
     fname = 'e1.arr'
     call write_tmp4

     tmp4(:,:,:) = wrk(1:nx,:,:,2) - wrk(1:nx,:,:,5)
     fname = 'e2.arr'
     call write_tmp4

     tmp4(:,:,:) = wrk(1:nx,:,:,3) - wrk(1:nx,:,:,6)
     fname = 'e3.arr'
     call write_tmp4



  end if
  return
end subroutine test_rhs_velocity
