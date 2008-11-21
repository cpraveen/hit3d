subroutine dealias_all

  use m_parameters
  use m_fields
  use x_fftw
  implicit none

  integer :: i,j,k

  ! truncating all the modes in which at least one component of the k-vector
  ! has magnitude that is larger than N/3

  do k = 1,nz
     do j = 1,ny
        do i = 1,nx+2
           if (sqrt(akx(i)**2+aky(k)**2+akz(j)**2).gt.dble(kmax) .or.  &
               akx(i).gt.nx/3  .or. aky(k).gt.nx/3 .or. akz(j).gt.nx/3) then
              fields(i,j,k,1:3+n_scalars) = zip
           end if
        end do
     end do
  end do

  return
end subroutine dealias_all