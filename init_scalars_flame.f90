
subroutine init_scalars

  use m_openmpi
  use m_io
  use m_parameters
  use m_fields
  use m_work
  use x_fftw

  !use LES

  implicit none

  integer :: k
  real*8  :: zloc, s, h0

  !h0 = 16.d0 * sqrt(tau_chem * REV)
  h0 = PI/8

  wrk(:,:,:,4) = dcmplx(0.d0,0.d0)

  ! creating 1/0 array of scalar
  do k = 1,nz
    zloc = dble(myid*nz + k-1) * dz
    s = 0.5*(tanh((zloc-PI*0.5)/h0) - tanh((zloc-PI*1.5)/h0));
    wrk(:,:,k,4) = dcmplx(s, 0.d0)
    !if (zloc.ge.PI/2.d0 .and. zloc.le.1.5d0*PI) then
    !  work(6,:,:,k) = dcmplx(1.d0,0.d0)
    !end if
  end do
  
  ! filtering the 1/0 array to make it smooth enough for numerical scheme
  !allocate(filter_g(nx,ny,nz))
  !call filter_fft_init(2,flame_thickness)
  !call filter_fft(6)
  !deallocate(filter_g)

  ! FFT of the filtered scalar
  call xFFT3d(1,4)

  ! putting it into the scalar array
  do k=1,n_scalars 
     fields(:,:,:,3+k) = wrk(:,:,:,4)
  end do

  return

end subroutine init_scalars
