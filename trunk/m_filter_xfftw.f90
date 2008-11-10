!======================================================================
!  Module that contains filtering rocedure for the part of the code that
!  addvects lagrangian particles
!
!  Time-stamp: <2008-11-06 11:19:23 (chumakov)>
!======================================================================
module m_filter_xfftw

  use m_parameters
  use x_fftw
  use m_work
  implicit none

  ! filter type (Gaussian by default)
  integer :: filter_type = 2

  ! filter size
  real*8 :: filter_size

  ! filter function in Fourier space
  real*8, allocatable :: filter_g(:,:,:)

!======================================================================
contains
!======================================================================

!======================================================================
!  filering subroutine for 2*pi^3-periodic domain
!  IN: ss(n,:,:,:)
!  OUT: ss(n:,;,:)
!======================================================================
  subroutine filter_xfftw(n)

    use m_io
    implicit none

    integer n

    integer :: ix, iy, iz

    real*8 :: a,b,c,d

!----------------------------------------------------------------------

    write(out,*) 'Filtering wrk #',n
    call flush(out)

!!$    call xFFT3d(1,n)

    ! multiplying wrk(:,:,:,n) by filter_g
    ! remember both are in the complex form
    do iy = 1, nz
       do iz = 1, nz_all
          do ix = 1, nx + 1, 2
             a = wrk(ix,iz,iy,n)
             b = wrk(ix+1,iz,iy,n)
             c = filter_g(ix,iz,iy)
             d = filter_g(ix+1,iz,iy)

             wrk(ix    , iz, iy, n) = a*c - b*d
             wrk(ix + 1, iz, iy, n) = b*c + a*d
          end do
       end do
    end do

!!$    ! perform inverse FFT
!!$    call xFFT3d(-1,n)

    return
  end subroutine filter_xfftw

!======================================================================
!  filering subroutine for 2*pi^3-periodic domain
!  IN: ss(n,:,:,:)
!  OUT: ss(n:,;,:)
!======================================================================
  subroutine filter_xfftw_fields(n)

    use m_io
    use m_fields
    implicit none

    integer n

    integer :: ix, iy, iz

    real*8 :: a,b,c,d

!----------------------------------------------------------------------

    write(out,*) 'Filtering fields #',n
    call flush(out)

!!$    call xFFT3d(1,n)

    ! multiplying wrk(:,:,:,n) by filter_g
    ! remember both are in the complex form
    do iy = 1, nz
       do iz = 1, nz_all
          do ix = 1, nx + 1, 2
             a = fields(ix,iz,iy,n)
             b = fields(ix+1,iz,iy,n)
             c = filter_g(ix,iz,iy)
             d = filter_g(ix+1,iz,iy)

             fields(ix    , iz, iy, n) = a*c - b*d
             fields(ix + 1, iz, iy, n) = b*c + a*d
          end do
       end do
    end do

!!$    ! perform inverse FFT
!!$    call xFFT3d(-1,n)

    return
  end subroutine filter_xfftw_fields

!======================================================================
!======================================================================
!======================================================================
!  subroutine that initializes the filter function filter_g
!  IN: filter_type - type of the filter
!      delta - filter's characteristic width
!
!  OUT: filter_g - normalized FFT of the filter function
!======================================================================
  subroutine filter_xfftw_init

    use m_parameters
    use m_io
    use m_work
    implicit none

    real*8  delta

    integer :: di,dj,dk,i,j,k, m
    integer :: idelta,sum1
    integer :: ii(8)
    real*8  :: rx,ry,rz
    real*8  :: const1,const2

    filter_size = particles_filter_size
    delta = filter_size

    if (delta .lt. dx*three) then
       write(out,*) "filter_xfftw_init: delta is too small:",delta,dx
       write(out,*) "Must be at least 3*dx = ",three*dx
       call my_exit(-1)
    end if

    write(out,*)'initializing the filter'

    ! looking for the available wrk array.  Later the filter will be
    ! created in this array and copied to filter_g array
!
    ! The logic is as follows.  At the end we need an array with the 
    ! Fourier transform of the filter function G; that array is filter_g.
    ! But our FFT routines work only on the "wrk" array, so to create 
    ! a full filter_g we need at least one full segment of the wrk array.
    ! By that we mean something like wrk(1:nx+2,1:ny,1:nz,n) for some n.
    ! The problem is that wrk array can be allocated with different size.
    ! So we get the dimensions of the currently allocated wrk, save in in the
    ! array ii(:), deallocate wrk, allocate it with the size that we need
    ! to define filter_g, use it to define filter_g, deallocate it and 
    ! allocate it again with the old sizes.

    ! note that since FFTW routines are tied to the size of wrk array,
    ! we need to deallocate/allocate the FFTW arrays and initialize them
    ! as well when we change the strufture of the wrk array

    ! getting the dimensions of currently allocated wrk array
    do i = 1,4
       ii(2*i-1) = LBOUND(wrk,i)
       ii(2*i  ) = UBOUND(wrk,i)
    end do
    write(out,"('Dimensions of wrk:',8i5)") ii(1:8)
    call flush(out)
    deallocate(wrk)
    allocate(wrk(1:nx+2,1:ny,1:nz,0:0),stat=ierr)
    if (ierr.ne.0) stop '*** m_filter_init: cannot allocate necessary wrk'
    wrk = zip

    ! re-initializing the FFTW arrays
    call x_fftw_allocate(-1)
    call x_fftw_allocate(1)
    call x_fftw_init

    ! number of the slice in the wrk array that will be used
    m = LBOUND(wrk,4)



!---------------------------------------------------------------
!   define the filtering function
!---------------------------------------------------------------
    case_filter_type: select case (filter_type)

!-------------------------------------------------------------------
!  tophat filter
!-------------------------------------------------------------------
    case (1)
       write(out,*) '-- tophat filter, delta =',delta
!       if(myid.eq.0) print*,'-- tophat filter, delta =',delta


       write(out,*) "Tophat filter is not working at the moment."
       write(out,*) "We're sorry for the inconvenience, stopping."
       call my_exit(-1)

       wrk(:,:,:,m) = zip

       idelta = delta / dx / 2
       ! normalization constant
       const1 = real((2*idelta+1)**3,8)
       const2 = 1.0d0 / const1
       do k = 1,nz
          dk = min(myid*nz+k-1,nz*numprocs-(myid*nz+k)+1)
          do j = 1,ny
             dj = min(j-1,ny-j+1)
             do i = 1,nx
                di = min(i-1,nx-i+1)
                wrk(i,j,k,m) = zip
                if (di.le.idelta .and. dj.le.idelta .and. dk.le.idelta) then
                   wrk(i,j,k,m) = const2
                end if
             end do
          end do
       end do
!-------------------------------------------------------------------
!  Gaussian filter
!-------------------------------------------------------------------
    case (2)
       write(out,*) '-- Gaussian filter, delta =',delta
       call flush(out)

       wrk(:,:,:,m) = zip

!       if(myid.eq.0) print*,'-- Gaussian filter, delta =',delta
       const1 = 6.0d0 / delta**2
       const2 = sqrt(const1/PI)**3 *dx**3
       do k = 1,nz
          dk = min(myid*nz+k-1,nz*numprocs-(myid*nz+k)+1)
          rz = dx * real(dk,8)
          do j = 1,ny
             dj = min(j-1,ny-j+1)
             ry = dx * real(dj,8)
             do i = 1,nx
                di = min(i-1,nx-i+1)
                rx = dx * real(di,8)
                rx = rx*rx+ry*ry+rz*rz
                wrk(i,j,k,m) = const2*exp(-const1*rx)

!                if (k==4) write(out,"(6i3,4e15.5)") i,j,k, di, dj, dk, rx, ry, rz, wrk(i,j,k,m)

             end do
          end do
       end do

       write(out,"(e15.6)") wrk(:,1,4,m)
       write(out,*) "------------------------"

!-------------------------------------------------------------------
!  linear filter
!-------------------------------------------------------------------
    case(3)
       write(out,*) '-- linear filter, delta =',delta
!       if(myid.eq.0) print*,'-- linear filter, delta =',delta
       const1 = 24.0d0 / (PI*delta**3)  *dx**3
       const2 = delta**2 / 4.0d0
       do k = 1,nz
          dk = min(myid*nz+k-1,nz*numprocs-(myid*nz+k)+1)
          rz = dx * real(dk,8)
          do j = 1,ny
             dj = min(j-1,ny-j+1)
             ry = dx * real(dj,8)
             do i = 1,nx
                di = min(i-1,nx-i+1)
                rx = dx * real(di,8)

                rx = rx*rx+ry*ry+rz*rz
                wrk(i,j,k,m) = zip
                if (rx.le.const2) wrk(i,j,k,m) = &
                     const1*(1.0d0-2.0d0*sqrt(rx)/delta)
             end do
          end do
       end do

    case default
       print *,'FILTER_FFT_INIT: wrong filter_type:',filter_type
       stop
    end select case_filter_type

    write(out,"(e15.6)") wrk(:,1,4,m)
    write(out,*) "------------------------"
    write(out,*) "------------------------"


    const1 = sum(wrk(1:nx,1:ny,1:nz,m))
    call MPI_REDUCE(const1,const2,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_TASK,mpi_err)
    if(myid.eq.0) write(out,*) '-- NORM OF G: ',const2

    write(out,"(e15.6)") wrk(:,1,4,m)
    write(out,*) "------------------------"
    write(out,*) "------------------------"

    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             tmp4(i,j,k) = wrk(i,j,k,m)
          end do
       end do
    end do
!!$    tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,m)
!!$    fname = 'blah'
!!$    call write_tmp4



!! outputting the filter
    write(out,"(e15.6)") wrk(:,1,4,m)


!---------------------------------------------------------------
!  compute FFT of g

    call xFFT3d(1,m)
!  print *,myid,' FILTER_FFT_INIT: FFT-1 done'

!  putting the FFT of the filter into filter_g.
    ! allocating the filter array
    if (.not.allocated(filter_g)) then
       allocate(filter_g(nx+2,ny,nz))
       filter_g = zip
       write(out,*) 'allocated filter_g'
       call flush(out)
    end if
    filter_g(:,:,:) = wrk(:,:,:,m)


    ! now we need to restore the size of the wrk array to what it was
    ! before the start of the subroutine
    deallocate(wrk)
    allocate(wrk(ii(1):ii(2),ii(3):ii(4),ii(5):ii(6),ii(7):ii(8)),stat=ierr)
    write(out,"('Re-allocated wrk(',i1,':',i4,',',i1,':',i4,',',i1,':',i4,',',i1,':',i4,')')") &
         ii(1:8)
    call flush(out)
    wrk = zip

    ! re-initializing the FFTW arrays
    call x_fftw_allocate(-1)
    call x_fftw_allocate(1)
    call x_fftw_init


!!$!  no need to normalize filter_g because our implementation of FFT
!!$!  normalizes the result anyways
!!$    const1 = one/real(nx**3,8)
!!$    do k=1,nz
!!$       do j=1,ny
!!$          do i=1,nx+2
!!$             filter_g(i,j,k) = wrk(i,j,k,m)*const1
!!$          end do
!!$       end do
!!$    end do

    return
  end subroutine filter_xfftw_init


end module m_filter_xfftw
