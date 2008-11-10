subroutine begin_restart

  use m_parameters
  use m_particles
  use m_fields
  use m_timing
  implicit none

  real*8, allocatable :: zhopa(:,:,:,:)
  integer :: i


  write(out,"(70('='))") 
  write(out,*) '                    RESTART '
  write(out,"(70('='))") 
  call flush(out)

  ITIME = ITMIN
  call get_file_ext

  ! reading the restart file
  if (task.eq.'hydro') call restart_read_parallel


!!$!================================================================================
!!$!  Checking the parallel read speed (reading 100 times)
!!$!  Currently the speedup factor from using parallel read is about 2.5
!!$!================================================================================
!!$  if(task.eq.'hydro') then
!!$     call m_timing_check
!!$     write(out,*) 'Start! ',cpu_min,cpu_sec
!!$     do i = 1,200
!!$        call restart_read! _parallel
!!$     end do
!!$     call m_timing_check
!!$     write(out,*) 'Finish!',cpu_min,cpu_sec
!!$
!!$  end if
!!$  call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!!$  stop 'checking the restart read speed'
!!$
!!$!================================================================================
!!$!  Checking the accuracy of the parallel read 
!!$!  (reading parallel first, then serial and comparing)
!!$!================================================================================
!!$  if (task.eq.'hydro') then
!!$
!!$     call restart_read_parallel
!!$
!!$     allocate(zhopa(nx+2,ny,nz,3+n_scalars), stat=ierr)
!!$     if (ierr.ne.0) stop 'cannot allocate zhopa'
!!$     zhopa = zip
!!$     zhopa = fields
!!$
!!$     call restart_read
!!$
!!$     print "(10e15.6)",maxval(abs(zhopa-fields)), fields(nx/2:nx/2+3,ny,nz,3), zhopa(nx/2:nx/2+3,ny,nz,3)
!!$  end if
!!$  call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!!$  stop 'checking the restart read'
!!$!================================================================================

  ! deciding whether we advance scalars or not
  if (time.gt.TSCALAR) then
     int_scalars = .true.
     write(out,"('Advancing ',i3,' scalars.')") n_scalars
  end if

  if (task.eq.'parts') call particles_init

  return
end subroutine begin_restart
