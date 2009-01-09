subroutine rhs_scalars

  use m_openmpi
  use m_io
  use m_parameters
  use m_fields
  use m_work
  use x_fftw
  use m_timing

  implicit none

  integer :: i, j, k, n, n1, n2, nv
  real*8  :: rtmp1, rtmp2, wnum2, r11, r12, r21, r22, r31, r32

  ! If we're not advancing scalars, put the real-space velocities 
  ! in wrk1...3 and return
  ! same is done if dealias=0, that is, we use 2/3 rule.
  if (.not.int_scalars) then
     ! converting velocities to the real space
     wrk(:,:,:,1:3) = fields(:,:,:,1:3)
     do n = 1,3
        call xFFT3d(-1,n)
     end do
     return
  end if

!--------------------------------------------------------------------------------
!  If dealias=0, performing the 2/3 rule dealiasing on scalars
!--------------------------------------------------------------------------------

  if (dealias.eq.0) then

     ! converting velocities to the real space
     wrk(:,:,:,1:3) = fields(:,:,:,1:3)
     do n = 1,3
        call xFFT3d(-1,n)
     end do

     ! Do each scalar one at a time.

     ! Trying to keep the velocities in wrk1:3 intact because they
     ! are needed later 
     do n = 1, n_scalars

        wrk(:,:,:,0) = fields(:,:,:,3+n)
        call xFFT3d(-1,0)

        ! Products of the scalar and velocities
        do i = 1,3
           wrk(:,:,:,n+2+i) = wrk(:,:,:,0) * wrk(:,:,:,i)
           call xFFT3d(1,n+2+i)
        end do

        ! Assembling the RHS in wrk(:,:,:,3+n)
        do k = 1,nz
           do j = 1,ny
              do i = 1,nx+1,2

                 ! If the dealiasing option is 2/3-rule (dealias=0) then we retain the modes
                 ! inside the cube described by $| k_i | \leq  k_{max}$, $i=1,2,3$.
                 ! The rest of the modes is purged

                 if (ialias(i,j,k) .gt. 0) then
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

                    wnum2 = akx(i)**2 + aky(k)**2 + akz(j)**2

                    wrk(i  ,j,k,3+n) =   rtmp1 - pe(n) * wnum2*fields(i  ,j,k,3+n)
                    wrk(i+1,j,k,3+n) = - rtmp2 - pe(n) * wnum2*fields(i+1,j,k,3+n)

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

  end if

!--------------------------------------------------------------------------------
!  If dealias=1, performing the phase shift and truncation on scalars
!  The main ideology is as follows.  First we evaluate the phase-shifted
!  quantities, then not phase shifted.  This is done in order to have more
!  working space: at the beginning we have all work arrays available, but at the
!  end we must put sin/cos factors in wrk0 and u.v.w in real space in wrk1...3.
!  They are going to be used again in rhs_velocity later.
!--------------------------------------------------------------------------------
  if (dealias.eq.1) then

     ! define the sin/cos factors that are used in phase shifting.
     ! computing sines and cosines for the phase shift of dx/2,dy/2,dz/2 
     ! and putting them into wrk0
     do k = 1,nz
        do j = 1,ny
           do i = 1,nx+1,2
              wrk(i  ,j,k,0) = cos(half*(akx(i  )+aky(k)+akz(j))*dx)
              wrk(i+1,j,k,0) = sin(half*(akx(i+1)+aky(k)+akz(j))*dx)
           end do
        end do
     end do

     ! Doing all scalars at the same time.  We can do this because the
     ! number of work arrays that are available to us is 3+n+2.  First
     ! three later will be taken by  velocities, and the last two are
     ! primary work arrays here.

     ! First, get phase shifted velocities and scalars
     do n = 1, 3+n_scalars
        ! phase-shifting the quantity
        do k = 1,nz
           do j = 1,ny
              do i = 1,nx+1,2
                 wrk(i  ,j,k,n) = fields(i  ,j,k,n) * wrk(i,j,k,0) - fields(i+1,j,k,n) * wrk(i+1,j,k,0)
                 wrk(i+1,j,k,n) = fields(i+1,j,k,n) * wrk(i,j,k,0) + fields(i  ,j,k,n) * wrk(i+1,j,k,0)
              end do
           end do
        end do
        ! transforming it to real space
        call xFFT3d(-1,n)
     end do

     ! now we have two vacant arrays: n_scalars+4 and n_scalars+5.  Work in them
     n1 = n_scalars+4; n2 = n_scalars+5;

     ! do one scalar at a time
     phase_shifted_rhs: do n = 4, 3+n_scalars

        ! getting all three products of phase-shifted scalar and phase-shifted velocities
        ! using three work arrays: n, n1 and n2
        wrk(:,:,:,n1) = wrk(:,:,:,n) * wrk(:,:,:,1)
        wrk(:,:,:,n2) = wrk(:,:,:,n) * wrk(:,:,:,2)
        wrk(:,:,:,n ) = wrk(:,:,:,n) * wrk(:,:,:,3)
        ! transforming them to Fourier space
        call xFFT3d(1,n1)  
        call xFFT3d(1,n2)  
        call xFFT3d(1,n )

        ! phase shifting them back using wrk0 and adding -0.5*(ik) to the RHS for the scalar
        do k = 1,nz
           do j = 1,ny
              do i = 1,nx+1,2

                 if (ialias(i,j,k) .gt. 1) then
                    wrk(i:i+1,j,k,n) = zip
                 else
                    ! (u+)*(phi+) phase shifted back
                    r11 = wrk(i  ,j,k,n1) * wrk(i,j,k,0) + wrk(i+1,j,k,n1) * wrk(i+1,j,k,0)
                    r12 = wrk(i+1,j,k,n1) * wrk(i,j,k,0) - wrk(i  ,j,k,n1) * wrk(i+1,j,k,0)
                    ! (v+)*(phi+) phase shifted back
                    r21 = wrk(i  ,j,k,n2) * wrk(i,j,k,0) + wrk(i+1,j,k,n2) * wrk(i+1,j,k,0)
                    r22 = wrk(i+1,j,k,n2) * wrk(i,j,k,0) - wrk(i  ,j,k,n2) * wrk(i+1,j,k,0)
                    ! (w+)*(phi+) phase shifted back
                    r31 = wrk(i  ,j,k,n ) * wrk(i,j,k,0) + wrk(i+1,j,k,n ) * wrk(i+1,j,k,0)
                    r32 = wrk(i+1,j,k,n ) * wrk(i,j,k,0) - wrk(i  ,j,k,n ) * wrk(i+1,j,k,0)
                    ! adding -0.5*(ik)*(the result) to the RHSs for the scalar
                    wrk(i  ,j,k,n) = + 0.5d0 * ( akx(i+1)*r12 + aky(k)*r22 + akz(j)*r32 )
                    wrk(i+1,j,k,n) = - 0.5d0 * ( akx(i  )*r11 + aky(k)*r21 + akz(j)*r31 )
                 end if

              end do
           end do
        end do

     end do phase_shifted_rhs

     ! at this moment wrk4...3+n_scalars contain the half of the convective term, which was
     ! obtained from the phase shifted quantities.  Now we need to add the other half of the
     ! convective term (the one that is obtained by multiplication of no-phase-shifted stuff)

     ! first get the velocities into the real space and put them in wrk1...3
     do n = 1,3
        wrk(:,:,:,n) = fields(:,:,:,n)
        call xFFT3d(-1,n)
     end do

     ! now do scalars one at a time since we don't have enough storage to do them 
     ! all at once

     not_phase_shifted_rhs: do n = 4, 3+n_scalars

        ! get the scalar into the real space and put it in the wrk(0)
        wrk(:,:,:,0) = fields(:,:,:,n)
        call xFFT3d(-1,0)

        ! <><><><><><><><><><><><> FOR NATA <><><><><><><><><><><><><><><><><><><><>
        ! This cycle is the right moment (probably) to add the reaction to the RHS
        ! since this is the only place where we have the scalar in real space.
        ! There are two choices:
        ! 1) do a fully dealiased reaction term. This requires adjustments in the 
        ! evaluation of phase-shifted stuff as well (you want to phase shift the
        ! scalar and compute s*(1-s) there, etc, and add 0.5 of it to RHS.  Then
        ! do the same here with not phase-shifted scalar.
        ! 2) compute the source term only here (with not phase-shifted scalar) but
        ! before you compute t*(1-t) in physical space you need to cut off all
        ! the modes outside of the small cube (all modes with ialias>0) from t.  This
        ! might require some thinking because later we need the scalar with the modes
        ! ialias<2, not ialias<1.  As a temporary measure, you can do it in wrk(n1)
        ! or wrk(n2) before they are used.
        ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


        ! calculate the product of the scalar with velocitiy components
        wrk(:,:,:,n1) = wrk(:,:,:,0) * wrk(:,:,:,1)
        wrk(:,:,:,n2) = wrk(:,:,:,0) * wrk(:,:,:,2)
        wrk(:,:,:, 0) = wrk(:,:,:,0) * wrk(:,:,:,3)
        ! transform it to the Fourier space
        call xFFT3d(1,n1)
        call xFFT3d(1,n2)
        call xFFT3d(1, 0)

        ! add the -0.5*(ik)*results to the RHS, along with the diffusion term
        do k = 1,nz
           do j = 1,ny
              do i = 1,nx+1,2

                 if (ialias(i,j,k) .lt. 2) then

                    rtmp1 =   akx(i+1)*wrk(i+1,j,k,n1) + aky(k)*wrk(i+1,j,k,n2) + akz(j)*wrk(i+1,j,k,0)
                    rtmp2 =   akx(i  )*wrk(i  ,j,k,n1) + aky(k)*wrk(i  ,j,k,n2) + akz(j)*wrk(i  ,j,k,0)

                    wnum2 = akx(i)**2 + aky(k)**2 + akz(j)**2
                    wrk(i  ,j,k,n) = wrk(i  ,j,k,n) + 0.5d0 * rtmp1 - pe(n-3) * wnum2*fields(i  ,j,k,n)
                    wrk(i+1,j,k,n) = wrk(i+1,j,k,n) - 0.5d0 * rtmp2 - pe(n-3) * wnum2*fields(i+1,j,k,n)

                 end if
              end do
           end do
        end do

     end do not_phase_shifted_rhs

  end if

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

           if ( abs(akx(i)).gt.akmax .or. &
                abs(aky(k)).gt.akmax .or. &
                abs(akz(j)).gt.akmax ) then

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
