module measurements
  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  use statistics
  implicit none

contains

  subroutine top_suscep(Q,Q2,susc_ave,susc_delta)
    real(dp), dimension(:), intent(in) :: Q, Q2
    real(dp), intent(out) :: susc_ave, susc_delta
    integer(i4) :: N,k,i
    real(dp) :: Qt,Q2t,jackk,Ntot
    real(dp), dimension(Mbins) :: Qmean,Q2mean,susc_avev
      N=size(Q)
      Ntot=real(N,dp)-real(N,dp)/real(Mbins,dp)
      call mean_0(Q,Qt)
      call mean_0(Q2,Q2t)
      susc_ave=Q2t-Qt**2
      Qmean=0._dp
      Q2mean=0._dp
      do i=1,Mbins
        do k=1,N
          if(k .le. (i-1)*N/Mbins) then
            Qmean(i)=Qmean(i)+Q(k)
            Q2mean(i)=Q2mean(i)+Q2(k)
          else if(k > i*N/Mbins) then
            Qmean(i)=Qmean(i)+Q(k)
            Q2mean(i)=Q2mean(i)+Q2(k)
          end if
        end do
        susc_avev(i)=(Q2mean(i)/Ntot) -(Qmean(i)/Ntot)**2
      end do
      do k=1,Mbins
        jackk=jackk+(susc_avev(k)-susc_ave )**2
      end do
      susc_delta=Sqrt(real(Mbins-1,dp)*jackk/real(Mbins,dp))
  end subroutine top_suscep

  subroutine initialize(corr1,corr2)
    real(dp), dimension(L,Nmsrs), intent(inout) :: corr1
    real(dp), dimension(L,Nmsrs), intent(inout) :: corr2
      corr1=0._dp
      corr2=0._dp
  end subroutine initialize

 subroutine correlation(Sx,Sy,k,corr1,corr2)
    real(dp), dimension(L), intent(in) :: Sx,Sy
    integer(i4), intent(in) :: k
    real(dp), dimension(L,Nmsrs), intent(inout) :: corr1
    real(dp), dimension(L,Nmsrs), intent(inout) :: corr2
    integer(i4) :: i1
    do i1=1,L
      corr1(i1,k)=Sx(i1)
      corr2(i1,k)=Sx(i1)*Sx(1)+Sy(i1)*Sy(1)
    end do
  end subroutine correlation

  subroutine correlation_function(corr1,corr2,CF,CFprom)
    real(dp), dimension(L,Nmsrs), intent(in) :: corr1
    real(dp), dimension(L,Nmsrs), intent(in) :: corr2
    real(dp), dimension(L), intent(out) :: CF,CFprom
    real(dp), dimension(L):: jackk
    real(dp) :: corr1m(L,Mbins)
    real(dp) :: corr2m(L,Mbins),CFm(L,Mbins)
    real(dp), dimension(L) :: corr1prom,corr1delta
    real(dp), dimension(L) :: corr2prom,corr2delta
    integer(i4) :: i1,i2,i3
    corr1prom=0._dp
    corr2prom=0._dp
    corr1delta=0._dp
    corr2delta=0._dp
    call mean_vector(corr1,corr1prom,corr1delta)
    call mean_vector(corr2,corr2prom,corr2delta)
    do i1=1,L
      CF(i1)=corr2prom(i1)-2._dp*corr1prom(i1)*corr1prom(1)
    end do
    jackk=0._dp
    do i1=1,Mbins
    corr1m(:,i1)=0._dp
    corr2m(:,i1)=0._dp
      do i2=1,Nmsrs
        if(i2 .le. (i1-1)*Nmsrs/Mbins) then
          corr1m(:,i1)=corr1m(:,i1)+corr1(:,i2)
          corr2m(:,i1)=corr2m(:,i1)+corr2(:,i2)
        else if(i2 > i1*Nmsrs/Mbins) then
          corr1m(:,i1)=corr1m(:,i1)+corr1(:,i2)
          corr2m(:,i1)=corr2m(:,i1)+corr2(:,i2)
        end if
      end do
    end do
    corr1m=corr1m/(real(Nmsrs,dp) -real(Nmsrs/Mbins,dp))
    corr2m=corr2m/(real(Nmsrs,dp) -real(Nmsrs/Mbins,dp))
    do i3=1,Mbins
      do i1=1,L
        CFm(i1,i3)=corr2m(i1,i3)-2._dp*corr1m(i1,i3)*corr1m(1,i3)
      end do
    end do
    do i1=1,Mbins
      jackk(:)=jackk(:)+(CF(:)-CFm(:,i1) )**2
    end do
    CFprom(:)=Sqrt(real(Mbins-1,dp)*jackk(:)/real(Mbins,dp))
  end subroutine correlation_function


end module measurements
