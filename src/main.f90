program main
  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use arrays
  use functions
  use statistics
  implicit none

  call thermalize(0.5_dp)
  !call vary_temp(0.1_dp,1._dp,30)
  !call test(0.5_dp)

contains

  subroutine thermalize(T)
    real(dp), intent(in) :: T
    real(dp), allocatable :: Sx(:),Sy(:)
    real(dp), allocatable :: ARtot(:)
    real(dp) :: AR,AR_ave,AR_delta
    integer(i4) :: i,k
    open(10, file = 'data/thermalization.dat', status = 'replace')
      allocate(Sx(L) )
      allocate(Sy(L))
      allocate(ARtot(Nmsrs))
      call hot_start(Sx,Sy)
      !call cold_start(Sx,Sy)
      k=0
      AR=0._dp
      do i=1,thermalization
        if(i==1 .or. mod(i,1)==0 ) then
          write(10,*) i, Hamilt(Sx,Sy)/real(L,dp)
        end if
        call Metropolis(T,Sx,Sy,AR)
        !call Cluster(T,Sx,Sy)
      end do

      do i=1,sweeps
        call Metropolis(T,Sx,Sy,AR)
        if(i>thermalization .and. mod(i,eachsweep)==0) then
          k=k+1
          ARtot(k)=AR
        end if
      end do
      call mean_scalar(ARtot,AR_ave,AR_delta)
      write(*,*) AR_ave,AR_delta
      deallocate(Sx,Sy,ARtot)
  end subroutine thermalize

  subroutine vary_temp(Ti,Tf,Nms)
    real(dp), intent(in) :: Ti,Tf
    integer(i4), intent(in) :: Nms
    real(dp), allocatable :: Sx(:),Sy(:)
    real(dp), allocatable :: ARtot(:),H(:),M2(:)
    real(dp) :: T,AR,AR_ave,AR_delta,H_ave,H_delta,M2_ave,M2_delta
    integer(i4) :: i,j,k
    open(10, file = 'data/accrate.dat', status = 'replace')
    open(20, file = 'data/energy.dat', status = 'replace')
    open(30, file = 'data/magnet2.dat', status = 'replace')
    do j=1,Nms
      dtheta=Pi*(0.1_dp+0.2_dp*real(j-1,dp)/real(Nms-1,dp))
      T=Ti+(Tf-Ti)*real(j-1,dp)/real(Nms-1,dp)
      allocate(Sx(L) )
      allocate(Sy(L))
      allocate(ARtot(Nmsrs))
      allocate(H(Nmsrs))
      allocate(M2(Nmsrs))
      !call hot_start(Sx,Sy)
      call cold_start(Sx,Sy)
      k=0
      !AR=0._dp
      do i=1,sweeps
        call Metropolis(T,Sx,Sy,AR)
        !call Cluster(T,Sx,Sy)
        if(i>thermalization .and. mod(i,eachsweep)==0) then
          k=k+1
          !write(*,*) k
          ARtot(k)=AR
          H(k)=Hamilt(Sx,Sy)
          M2(k)=Magnet2(Sx,Sy)
        end if
      end do
      call mean_scalar(ARtot,AR_ave,AR_delta)
      call mean_scalar(H,H_ave,H_delta)
      call mean_scalar(M2,M2_ave,M2_delta)
      write(10,*) T,AR_ave,AR_delta
      write(20,*) T,H_ave/real(L,dp),H_delta/real(L,dp)
      write(30,*) T,M2_ave/real(L,dp),M2_delta/real(L,dp)
      deallocate(Sx,Sy,ARtot,H,M2)
    end do
    close(10)
    close(20)
    close(30)
  end subroutine vary_temp

  subroutine test(T)
    real(dp), intent(in) :: T
    real(dp), allocatable :: Sx(:),Sy(:)
    integer(i4) :: i
      allocate(Sx(L) )
      allocate(Sy(L))
      call hot_start(Sx,Sy)
      !call cold_start(Sx,Sy)
      do i=1,10
        call Clustert(T,Sx,Sy)
        write(*,*) "Is it 1?", i, Sx(1)**2+Sy(1)**2
      end do
  end subroutine test

  subroutine test2(T)
    real(dp), intent(in) :: T
    real(dp), allocatable :: Sx(:),Sy(:)
    integer(i4) :: i
    open(10, file = 'data/thermalization.dat', status = 'replace')
      allocate(Sx(L) )
      allocate(Sy(L))
      call hot_start(Sx,Sy)
      !call cold_start(Sx,Sy)
      do i=1,thermalization
        call Cluster(T,Sx,Sy)
        if(i==1 .or. mod(i,eachsweep)==0 ) then
          write(10,*) i, Hamilt(Sx,Sy)/real(L,dp)
        end if
      end do
    close(10)
  end subroutine test2

end program main
