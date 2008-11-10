!================================================================================
! Module contains interface to OpenMPI
!
! Time-stamp: <2008-09-09 11:17:57 (chumakov)>
!================================================================================


module m_openmpi
!================================================================================
  implicit none
  include 'mpif.h'

  ! --- MPI variables
  logical :: iammaster
  integer(kind=MPI_INTEGER_KIND) :: myid_world, numprocs_world
  integer(kind=MPI_INTEGER_KIND) :: numprocs_hydro, numprocs_stats, numprocs_parts
  integer(kind=MPI_INTEGER_KIND) :: myid, numprocs, master, mpi_err, mpi_info
  integer(kind=MPI_INTEGER_KIND) :: id_to, id_from, tag, count
  integer(kind=MPI_INTEGER_KIND) :: id_root_hydro, id_root_stats, id_root_parts

  ! communicator for separate tasks
  integer(kind=MPI_INTEGER_KIND) :: MPI_COMM_TASK
  ! exclusive communicator for root processes of tasks
  integer(kind=MPI_INTEGER_KIND) :: MPI_COMM_ROOTS

  integer (kind=MPI_INTEGER_KIND) :: sendtag, recvtag
  integer (kind=MPI_INTEGER_KIND) :: request, request1, request2, request3, mpi_request
  integer (kind=MPI_INTEGER_KIND) :: id_l, id_r
  integer (kind=mpi_INTEGER_KIND) :: mpi_status(MPI_STATUS_SIZE)

  integer(kind=MPI_INTEGER_KIND) :: color, key

  character*5 :: task
  character*10 :: run_name_local

!================================================================================
contains
!================================================================================

  subroutine m_openmpi_init

    implicit none

    integer (kind=mpi_INTEGER_KIND) :: n
    integer*4 :: np_local
    integer :: i

    ! first getting the run name (it's local, not  global run_name)
    call openmpi_get_run_name


    ! initializing MPI environment
    call MPI_INIT(mpi_err)
    call MPI_Comm_size(MPI_COMM_WORLD,numprocs_world,mpi_err)
    call MPI_Comm_rank(MPI_COMM_WORLD,myid_world,mpi_err)

!--------------------------------------------------------------------------------
!  splitting the communicator into several parts (hydro, stat, particles etc)
!--------------------------------------------------------------------------------

    ! first finding out if there are any particles involved.
    ! if there are no particles, then we split the processors in two parts:
    ! hydro and stats.  If there are some particles, we split the processors
    ! in three parts: "hydro", "stats" and "parts".  The variables that 
    ! determines which part the process belongs to is "task".

    ! first see, how many particles are there
    if (myid_world.eq.0) then
       ! opening the inupt file
       open(99,file=run_name_local//'.in')
       ! skipping the first 35 lines
       do i = 1,35
          read(99,*)
       end do
       ! reading the number of particles
       read(99,*) np_local
       close(99)
    end if

    ! broadcasting the number of particles to all processors
    count = 1
    call MPI_BCAST(np_local,count,MPI_INTEGER4,0,MPI_COMM_WORLD,mpi_err)

    ! now splitting the processors in tasks: hydro, stats and parts
    ! the curren logic is this:
    ! - In case if there are no particles, the split between the hydro and
    ! the stats part is 2/3 and 1/3.  This way the total number of processors
    ! needs to be 3*2^n
    ! - In case with the particles in the flow, the split is 1/2, 1/4 and 1/4

    ! first the case when we do not have particles
    if (np_local.eq.0) then

       ! if the numprocs_total is divisible by 3, assign 2/3 of it to hydro
       ! and the rest to stats
       if (int(numprocs_world/3)*3 .eq. numprocs_world) then
          numprocs_hydro = numprocs_world * 2/3
          numprocs_stats = numprocs_world - numprocs_hydro
          numprocs_parts = 0
       else if (2**floor(log(real(numprocs_world))/log(2.d0)) .eq. numprocs_world) then
          print*, 'numprocs_world is 2^n, allocating half for hydro: ',numprocs_world
          numprocs_hydro = numprocs_world / 2
          numprocs_stats = numprocs_world - numprocs_hydro
          numprocs_parts = 0
       else
          ! if the # of processors N is not 2^n and not divisible by 3, then just take
          ! the biggest 2^k < N and make these hydro, the rest - stat.
          numprocs_hydro = 2**floor(log(real(numprocs_world))/log(2.d0))
          numprocs_stats = numprocs_world - numprocs_hydro
          numprocs_parts = 0
       end if

       id_root_hydro = 0
       id_root_stats = numprocs_hydro
       id_root_parts = 0

    else
       numprocs_hydro = numprocs_world / 2
       numprocs_stats = numprocs_world / 4
       numprocs_parts = numprocs_world / 4

       id_root_hydro = 0
       id_root_stats = numprocs_hydro
       id_root_parts = numprocs_hydro + numprocs_stats

    end if

    ! splitting the communicator into several parts
    ! (currently two)
    ! 1. hydro
    ! 2. stats
    ! 3. parts

    if (myid_world.lt.numprocs_hydro) then
       task = 'hydro'
       color = 0
       myid = myid_world
    elseif (myid_world.ge.numprocs_hydro .and. myid_world .lt. numprocs_hydro+numprocs_stats) then
       task = 'stats'
       color = 1
       myid = myid_world - numprocs_hydro
    else
       task = 'parts'
       color = 2
       myid = myid_world - numprocs_hydro - numprocs_stats
    end if
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,myid,MPI_COMM_TASK,mpi_err)
    call MPI_COMM_SIZE(MPI_COMM_TASK,numprocs,mpi_err)
    call MPI_COMM_RANK(MPI_COMM_TASK,myid,mpi_err)


    ! seeing if this is the master process
    master = 0
    iammaster = .false.
    if (myid.eq.master) iammaster=.true.


!!$    ! The following is put on hold because it looks like a crazy idea

!!$    ! now creating separate exclusive communicator for the master nodes only 
!!$    ! the name of the new communicator is MPI_COMM_ROOTS
!!$    ! if we want quickly broadcast something, then we can use two BCAST calls
!!$    color = 1
!!$    if (iammaster) color = 0
!!$    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,myid_world,MPI_COMM_ROOTS,mpi_err)



    return
  end subroutine m_openmpi_init

!================================================================================

  subroutine m_openmpi_exit

    call MPI_COMM_FREE(MPI_COMM_TASK,mpi_err)
    call MPI_FINALIZE(mpi_err)

    return
  end subroutine m_openmpi_exit


!================================================================================

  subroutine openmpi_get_run_name
    implicit none

    character*80 :: tmp_str
    integer      :: iargc

    ! reading the run_name from the command line
    if(iargc().eq.0) then
       call getarg(0,tmp_str)
       print*, 'Format: ',trim(tmp_str),' <run name>'
       stop
    end if
    call getarg(1,run_name_local)
    if(len_trim(run_name_local).ne.10) then
       print *, 'Run name: "',run_name_local,'"'
       print *, '          "1234567890"'
       print *, 'Length of run name is less than 10, sorry.'
       stop
    end if

  end subroutine openmpi_get_run_name

!================================================================================

end module m_openmpi
