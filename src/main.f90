program main
  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use arrays
  use functions
  use statistics
  implicit none

  !call thermalize(0.2_dp)
  !call vary_temp(0.1_dp,2._dp,30)
  call fixed_temp(0.5_dp)
  !call test(0.1_dp)

contains

  subroutine thermalize(T)
    real(dp), intent(in) :: T
    real(dp), allocatable :: Sx(:),Sy(:)
    real(dp), allocatable :: ARtot(:)
    real(dp) :: AR,AR_ave,AR_delta
    integer(i4) :: i,k
    open(10, file = 'data/ThE.dat', status = 'replace')
    open(20, file = 'data/ThQ.dat', status = 'replace')
      allocate(Sx(L) )
      allocate(Sy(L))
      allocate(ARtot(Nmsrs))
      call hot_start(Sx,Sy)
      !call cold_start(Sx,Sy)
      k=0
      AR=0._dp
      do i=1,2*thermalization
        if(i==1 .or. mod(i,1)==0 ) then
          write(10,*) i, Hamilt(Sx,Sy)/real(L,dp)
          write(20,*) i, top_charge(Sx,Sy)
        end if
        !call Metropolis(T,Sx,Sy,AR)
        call Cluster(T,Sx,Sy,AR)
      end do

      do i=1,sweeps
        call Metropolis(T,Sx,Sy,AR)
        if(i>thermalization .and. mod(i,eachsweep)==0) then
          k=k+1
          ARtot(k)=AR
        end if
      end do
      call mean_scalar(ARtot,AR_ave,AR_delta)
      write(*,*) "Acc. Rate=",AR_ave,AR_delta
      close(10)
      close(20)
      deallocate(Sx,Sy,ARtot)
  end subroutine thermalize

  subroutine vary_temp(Ti,Tf,Nms)
    real(dp), intent(in) :: Ti,Tf
    integer(i4), intent(in) :: Nms
    real(dp), allocatable :: Sx(:),Sy(:)
    real(dp), allocatable :: ARtot(:),H(:),M2(:),Q(:),Q2(:)
    real(dp) :: T,AR,AR_ave,AR_delta,H_ave,H_delta,M2_ave,M2_delta
    real(dp) :: Q_ave,Q_delta,Q2_ave,Q2_delta
    integer(i4) :: i,j,k
    open(10, file = 'data/accrate.dat', status = 'replace')
    open(20, file = 'data/energy.dat', status = 'replace')
    open(30, file = 'data/magnet2.dat', status = 'replace')
    open(40, file = 'data/topcharge.dat', status = 'replace')
    open(50, file = 'data/topcharge2.dat', status = 'replace')
    do j=1,Nms
      dtheta=Pi*(0.1_dp+0.5_dp*real(j-1,dp)/real(Nms-1,dp))
      T=Ti+(Tf-Ti)*real(j-1,dp)/real(Nms-1,dp)
      allocate(Sx(L) )
      allocate(Sy(L))
      allocate(ARtot(Nmsrs))
      allocate(H(Nmsrs))
      allocate(M2(Nmsrs))
      allocate(Q(Nmsrs))
      allocate(Q2(Nmsrs))
      call hot_start(Sx,Sy)
      !call cold_start(Sx,Sy)
      k=0
      !AR=0._dp
      do i=1,sweeps
        !call Metropolis(T,Sx,Sy,AR)
        call Cluster(T,Sx,Sy,AR)
        if(i>thermalization .and. mod(i,eachsweep)==0) then
          k=k+1
          !write(*,*) k
          ARtot(k)=AR
          H(k)=Hamilt(Sx,Sy)
          M2(k)=Magnet2(Sx,Sy)
          Q(k)=top_charge(Sx,Sy)
          Q2(k)=Q(k)**2
        end if
      end do
      call mean_scalar(ARtot,AR_ave,AR_delta)
      call mean_scalar(H,H_ave,H_delta)
      call mean_scalar(M2,M2_ave,M2_delta)
      call mean_scalar(Q,Q_ave,Q_delta)
      call mean_scalar(Q2,Q2_ave,Q2_delta)
      write(10,*) T,AR_ave,AR_delta
      write(20,*) T,H_ave/real(L,dp),H_delta/real(L,dp)
      write(30,*) T,M2_ave/real(L,dp),M2_delta/real(L,dp)
      write(40,*) T,Q_ave,Q_delta
      write(50,*) T,Q2_ave,Q2_delta
      deallocate(Sx,Sy,ARtot,H,M2,Q,Q2)
    end do
    close(10)
    close(20)
    close(30)
    close(40)
    close(50)
  end subroutine vary_temp

  subroutine test(T)
    real(dp), intent(in) :: T
    real(dp), allocatable :: Sx(:),Sy(:)
    real(dp) :: x
    integer(i4) :: i
      allocate(Sx(L) )
      allocate(Sy(L))
      call hot_start(Sx,Sy)
      !call cold_start(Sx,Sy)
      do i=1,10
        call Cluster(T,Sx,Sy,x)
        write(*,*) "Is it 1?", i, Sx(1)**2+Sy(1)**2
      end do
  end subroutine test

  subroutine fixed_temp(T)
    real(dp), intent(in) :: T
    real(dp), allocatable :: Sx(:),Sy(:)
    real(dp), allocatable :: result2(:)
    integer(i4), allocatable :: result1(:)
    real(dp) :: x,Q,norm
    integer(i4) :: i
    open(10, file = 'data/histogram.dat', status = 'replace')
    allocate(Sx(L) )
    allocate(Sy(L))
    allocate(result1(bins))
    allocate(result2(bins))
    do i=1,bins
      result2(i) =minn+binwidth/2._dp+real(i-1,dp)*binwidth
    end do
    result1=0

    call hot_start(Sx,Sy)
    !call cold_start(Sx,Sy)
    do i=1,sweeps
      call Cluster(T,Sx,Sy,x)
      if(i>thermalization .and. mod(i,eachsweep)==0 ) then
        Q=top_charge(Sx,Sy)
        call histogram(Q,result1,result2)
      end if
    end do

    norm=0._dp
    do i=1,bins
      norm=norm+result1(i)
    end do
    !norm=norm*4._dp/101._dp

    do i=1,bins
      write(10,*) result2(i), result1(i)/norm, sqrt( real(result1(i),dp) )/norm
    end do

    close(10)
  end subroutine fixed_temp

end program main
