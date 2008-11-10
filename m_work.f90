!================================================================================
! M_WORK - module that contains working arrays wrk1....wrk15.
!
! Time-stamp: <2008-11-06 11:14:39 (chumakov)>
!================================================================================
module m_work

!  use m_openmpi
  use m_parameters
  use m_io


  implicit none

  ! --- work arrays
  real*4,allocatable :: tmp4(:,:,:)

  real*8, allocatable :: wrk(:,:,:,:)

  real*8, allocatable :: wrkp(:,:)

  real*8, allocatable :: rhs_old(:,:,:,:)

contains

!================================================================================
  subroutine m_work_init

    use m_parameters
    implicit none

    ! allocating work arrays
    if (task.eq.'hydro')  then
       call m_work_allocate(max(6,5+n_scalars))
    elseif (task.eq.'stats')  then
       call m_work_allocate(6)
    elseif (task.eq.'parts') then
       ! required recources are different for the "parts" part of the code
       ! we need several layers of velocities for velocity interpolation
       select case (particles_tracking_scheme)
       case (0)
          allocate(wrk(1:nx+2,1:ny,1:3,0:0),stat=ierr)
          write(out,*) "Allocated wrk(1:nx+2,1:ny,1:3,0:0)"
          wrk = zip
       case (1)
          allocate(wrk(1:nx+2,1:ny,1:3,1:3),stat=ierr)
          write(out,*) "Allocated wrk(1:nx+2,1:ny,1:3,1:3)"
          wrk = zip
       case default
          stop 'wrong particles_tracking_scheme'
       end select
       allocate(tmp4(nx,ny,nz), stat=ierr)
       write(out,*) "Allocated tmp4."
       call flush(out)

    else
       print *,'TASK variable is set in such a way that I dont know how to allocate work arrays'
       print *,'task = ',task
    end if


    return
  end subroutine m_work_init

!================================================================================

  ! --- allocating and zeroing out the prescribed number of arrays
  subroutine m_work_allocate(number)

    implicit none
    integer :: number
    integer :: i,ierr

    ierr = 0

    write(out,"('Allocating work: ',i3)") number
    call flush(out)

    ! array that is needed for output (nx,ny,nz)
    allocate(tmp4(nx,ny,nz),stat=i);   ierr = ierr + i;

    ! main working array, needed for FFT etc, so (nx+2,ny,nz)
    allocate(wrk(nx+2,ny,nz,0:number),stat=i);  ierr = ierr + i

    ! array for the spare RHS for Adams-Bashforth time-stepping scheme methods
    if (task.eq.'hydro') then
       allocate(rhs_old(nx+2,ny,nz,3+n_scalars),stat=i); ierr = ierr + i
    end if

    if (ierr.ne.0) then
       print *,'*** WORK_INIT: error in allocation, stopping. ',ierr
       print *,'*** task = ',task
       print *,'*** myid = ',myid
       call my_exit(-1)
       stop
    end if

    write(out,"('allocated wrk(nx+2,ny,nz,0:',i3,')')") number
    write(out,"('allocated rhs_old(nx+2,ny,nz,1:',i3,')')") 3+n_scalars
    call flush(out)



    tmp4 = 0.0
    wrk = 0.0d0
    rhs_old = 0.d0


    return
  end subroutine m_work_allocate


  subroutine m_work_exit
    implicit none

!    write(out,*) 'deallocaing tmp4';  call flush(out)
    if (allocated(tmp4)) deallocate(tmp4)   
!    write(out,*) 'deallocaing wrk';  call flush(out)
    if (allocated(wrk)) deallocate(wrk)   
!    write(out,*) 'deallocaing wrkp';  call flush(out)
    if (allocated(wrkp)) deallocate(wrkp)
    write(out,*) 'deallocated wrk arrays';  call flush(out)

    return
  end subroutine m_work_exit

end module m_work
