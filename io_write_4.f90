subroutine io_write_4

  ! Writing out the velocities and scalars in X-space
  ! to the real*4 file

  use m_parameters
  use m_fields
  use m_work
  use x_fftw

  implicit none

  integer :: n

  ! velocities
  wrk(:,:,:,1) = fields(:,:,:,1)
  call xFFT3d(-1,1)
  fname = 'u.'//file_ext
  tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,1)
!!$  write(out,*) "Writing u with all procs"
!!$  call flush(out)
  call write_tmp4_all

!!$  fname = 'u1.'//file_ext
!!$  write(out,*) "Writing u with one proc"
!!$  call flush(out)
!!$  call write_tmp4
!!$  write(out,*) "wrote."
!!$  call flush(out)

  wrk(:,:,:,1) = fields(:,:,:,2)
  call xFFT3d(-1,1)
  fname = 'v.'//file_ext
  tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,1)
  call write_tmp4_all

  wrk(:,:,:,1) = fields(:,:,:,3)
  call xFFT3d(-1,1)
  fname = 'w.'//file_ext
  tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,1)
  call write_tmp4_all

  ! scalars
  if (int_scalars) then
     do n = 1,n_scalars
        wrk(:,:,:,1) = fields(:,:,:,3+n)
        call xFFT3d(-1,1)
        write(fname,"('sc',i2.2,'.',a6)") n,file_ext
        tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,1)
        call write_tmp4_all

     end do
  end if

  return
end subroutine io_write_4
