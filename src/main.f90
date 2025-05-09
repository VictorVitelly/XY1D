program main
  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use arrays
  use functions
  use statistics
  use measurements
  implicit none

  call cpu_time(starting)
  !call thermalize(0.2_dp)
  !call vary_temp(0.1_dp,2._dp,20)
  !call fixed_temp(0.5_dp)
  !call test(0.1_dp)
  call correlate(0.1_dp,1.0_dp,10)

  call cpu_time(ending)
  write(*,*) "Elapsed time: ", (ending-starting), " s"

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
    real(dp), allocatable :: ARtot(:),H(:),M2(:),Q(:),Q2(:),Q4(:)
    real(dp) :: T,AR,AR_ave,AR_delta,H_ave,H_delta,M2_ave,M2_delta
    real(dp) :: Q_ave,Q_delta,Q2_ave,Q2_delta,susc_ave,susc_delta
    real(dp) :: Q4_ave,Q4_delta,c4,c4_delta
    integer(i4) :: i,j,k
    open(10, file = 'data/accrate.dat', status = 'replace')
    open(20, file = 'data/energy.dat', status = 'replace')
    open(30, file = 'data/magnet2.dat', status = 'replace')
    open(40, file = 'data/topcharge.dat', status = 'replace')
    open(50, file = 'data/topcharge2.dat', status = 'replace')
    open(60, file = 'data/c4.dat', status = 'replace')
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
      allocate(Q4(Nmsrs))
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
          Q4(k)=Q(k)**4
        end if
      end do
      call mean_scalar(ARtot,AR_ave,AR_delta)
      call mean_scalar(H,H_ave,H_delta)
      call mean_scalar(M2,M2_ave,M2_delta)
      call mean_scalar(Q,Q_ave,Q_delta)
      call mean_scalar(Q2,Q2_ave,Q2_delta)
      call mean_scalar(Q4,Q4_ave,Q4_delta)
      call top_suscep(Q,Q2,susc_ave,susc_delta)
      write(10,*) T,AR_ave,AR_delta
      write(20,*) T,H_ave/real(L,dp),H_delta/real(L,dp)
      write(30,*) T,M2_ave/real(L,dp),M2_delta/real(L,dp)
      write(40,*) T,Q_ave,Q_delta
      write(50,*) T,Q2_ave/real(L,dp),Q2_delta/real(L,dp),susc_ave/real(L,dp),susc_delta/real(L,dp)
      c4=(3._dp*(Q2_ave**2)-Q4_ave)/real(L,dp)
      c4_delta=sqrt((6._dp*(Q2_ave*Q2_delta)**2+Q4_delta**2)/real(L,dp))
      write(60,*) T,c4, c4_delta
      deallocate(Sx,Sy,ARtot,H,M2,Q,Q2,Q4)
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
    do i=1,bins
      write(10,*) result2(i), result1(i)/norm, sqrt( real(result1(i),dp) )/norm
    end do
    close(10)
  end subroutine fixed_temp

  subroutine correlate(T0,Tf,NTs)
  real(dp), intent(in) :: T0,Tf
  integer(i4), intent(in) :: NTs
  integer(i4) :: i,k,j,k2
  real(dp) :: Sx(L),Sy(L)
  real(dp) :: corr1(L,Nmsrs),corr2(L,Nmsrs),CF(L),CFprom(L)
  real(dp), allocatable :: results(:,:),deltaresults(:,:),Q2(:)
  real(dp) :: T,x,Q2_ave,Q2_delta
  open(60, file = 'data/corrfunc.dat', status = 'replace')
  open(50, file = 'data/topsuscep.dat', status = 'replace')
    allocate(results(L,NTs) )
    allocate(deltaresults(L,NTs) )
    k2=0
    do j=1,NTs
      allocate(Q2(Nmsrs))
      T=T0+(Tf-T0)*real(j-1,dp)/real(NTs-1,dp)
      write(*,*) T
      k2=k2+1
      call hot_start(Sx,Sy)
      call initialize(corr1,corr2)
      k=0
      do i=1,sweeps
        call Cluster(T,Sx,Sy,x)
        if(i>thermalization .and. mod(i,eachsweep)==0) then
          k=k+1
          Q2(k)=top_charge(Sx,Sy)**2
          call correlation(Sx,Sy,k,corr1,corr2)
        end if
      end do
      call mean_scalar(Q2,Q2_ave,Q2_delta)
      write(50,*) T,Q2_ave/real(L,dp),Q2_delta/real(L,dp)

      call correlation_function(corr1,corr2,CF,CFprom)
      do i=1,L
        results(i,k2)=CF(iv(i))
        deltaresults(i,k2)=CFprom(iv(i))
      end do
      deallocate(Q2)
    end do
    do i=1,L
      write(60,*) abs(i-1), results(i,:), deltaresults(i,:)
    end do
    close(50)
    close(60)
    deallocate(results,deltaresults)
  end subroutine correlate

end program main
