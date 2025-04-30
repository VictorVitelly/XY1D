module statistics
  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  implicit none

contains

  subroutine random_real(x,bound)
    real(dp),intent(out) :: x
    real(dp), intent(in) :: bound
    real(dp) :: y
    call random_number(y)
    x = 2._dp*bound*y -bound
  end subroutine random_real

  subroutine cold_start(Sx,Sy)
    real(dp), dimension(:), intent(out) :: Sx,Sy
    Sx=1.0_dp
    Sy=0.0_dp
  end subroutine cold_start

  subroutine hot_start(Sx,Sy)
    real(dp), dimension(:), intent(out) :: Sx,Sy
    integer(i4) :: i,Narr
    real(dp) :: theta,x,y
    Narr=size(Sx)
    do i=1,Narr
      call random_real(theta,Pi)
      x=cos(theta)
      y=sin(theta)
      Sx(i)=x/(x**2+y**2)
      Sy(i)=y/(x**2+y**2)
    end do
  end subroutine hot_start

  subroutine Metropolis(T,Sx,Sy,AR)
    real(dp), intent(in) :: T
    real(dp), dimension(:), intent(inout) :: Sx,Sy
    real(dp),intent(out) :: AR
    integer(i4) :: i,Narr,AR_N
    real(dp) :: r,theta,Sx2,Sy2,DE,p,a,b
    Narr=size(Sx)
    AR=0._dp
    AR_N=0
    do i=1,Narr
      call random_real(theta,dtheta)
      a=cos(theta)*Sx(i)-sin(theta)*Sy(i)
      b=sin(theta)*Sx(i)+cos(theta)*Sy(i)
      Sx2=a/(a**2+b**2)
      Sy2=b/(a**2+b**2)
      DE=DeltaE(T,Sx,Sy,i,Sx2,Sy2)
      if(DE.le.0._dp) then
        Sx(i)=Sx2
        Sy(i)=Sy2
        AR=AR+1._dp
        AR_N=AR_N+1
      else
        call random_number(r)
        p=Exp(-DE)
        AR=AR+p
        AR_N=AR_N+1
        if(r<p) then
          Sx(i)=Sx2
          Sy(i)=Sy2
        end if
      end if
    end do
    AR=AR/real(AR_N,dp)
  end subroutine Metropolis

  subroutine Cluster(T,Sx,Sy)
    real(dp), intent(in) :: T
    real(dp), dimension(:), intent(inout) :: Sx,Sy
    real(dp) rx,ry,theta,r,dh,p,p2
    real(dp), allocatable :: Sx2(:),Sy2(:)
    integer(i4) :: i,j,Narr,k1,k2
    integer(i4), allocatable :: bond(:),clusterL(:),clusterR(:),clusterRaux(:)
    Narr=size(Sx)
    allocate(bond(Narr))
    allocate(Sx2(Narr))
    allocate(Sy2(Narr))
    allocate(clusterL(Narr))
    allocate(clusterR(Narr))
    allocate(clusterRaux(Narr))
    call random_real(theta,Pi)
    rx=cos(theta)
    ry=sin(theta)
    bond=0
    do i=1,Narr
      Sx2(i)=Sx(i)-2._dp*rx*(Sx(i)*rx+Sy(i)*ry )
      Sy2(i)=Sy(i)-2._dp*ry*(Sx(i)*rx+Sy(i)*ry )
      !Sx2(i)=Sx2(i)/(Sx2(i)**2+Sy2(i)**2)
      !Sy2(i)=Sy2(i)/(Sx2(i)**2+Sy2(i)**2)
    end do
    do i=1,Narr
      dh=Deltah(T,Sx,Sy,i,Sx2(iv(i+1)),Sy2(iv(i+1)))
      if(dh > 0) then
        call random_number(r)
        p=1._dp-Exp(-dh/T)
        if(r<p) then
          bond(i)=1
        end if
      end if
    end do
    k1=0
    k2=0
    do i=1,Narr
      if(bond(iv(i-1))==0 .and. bond(i)==1 ) then
        k1=k1+1
        clusterL(k1)=i
      end if
      if(bond(iv(i+1))==0 .and. bond(i)==1 ) then
        k2=k2+1
        clusterR(k2)=i
      end if
    end do
    if(k1.ne.k2) then
      write(*,*) "WARNING!"
    end if
    if(clusterL(1)>clusterR(1) ) then
      do i=1,k2-1
        clusterRaux(i)=clusterR(iv(i+1) )
      end do
      clusterRaux(k2)=clusterR(1)
      clusterR(:)=0._dp
      clusterR(:)=clusterRaux(:)
    end if
    do i=1,k1
      call random_number(p2)
      if(p2<0.5_dp) then
        do j=clusterL(i),clusterR(i)
          Sx(j)=Sx2(j)/(Sx2(j)**2 +Sy2(j)**2 )
          Sy(j)=Sy2(j)/(Sx2(j)**2 +Sy2(j)**2 )
        end do
      end if
    end do
    deallocate(bond,Sx2,Sy2,clusterL,clusterR,clusterRaux)
  end subroutine Cluster


  subroutine Clustert(T,Sx,Sy)
    real(dp), intent(in) :: T
    real(dp), dimension(:), intent(inout) :: Sx,Sy
    real(dp) rx,ry,theta,r,dh,p
    real(dp), allocatable :: Sx2(:),Sy2(:)
    integer(i4) :: i,j,Narr,k1,k2
    integer(i4), allocatable :: bond(:),clusterL(:),clusterR(:),clusterRaux(:)
    Narr=size(Sx)
    allocate(bond(Narr))
    allocate(Sx2(Narr))
    allocate(Sy2(Narr))
    allocate(clusterL(Narr))
    allocate(clusterR(Narr))
    allocate(clusterRaux(Narr))
    call random_real(theta,Pi)
    rx=cos(theta)
    ry=sin(theta)
    bond=0
    do i=1,Narr
      Sx2(i)=Sx(i)-2._dp*rx*(Sx(i)*rx+Sy(i)*ry )
      Sy2(i)=Sy(i)-2._dp*ry*(Sx(i)*rx+Sy(i)*ry )
    end do
    do i=1,Narr
      dh=Deltah(T,Sx,Sy,i,Sx2(iv(i+1)),Sy2(iv(i+1)))
      if(dh > 0) then
        call random_number(r)
        p=1._dp-Exp(-dh/T)
        if(r<p) then
          bond(i)=1
        end if
      end if
    end do
    write(*,*) "bound=", bond
    k1=0
    k2=0
    do i=1,Narr
      if(bond(i)==0 ) then
        k1=k1+1
        clusterL(k1)=i
        clusterR(k1)=i
      else if(bond(i)==1) then
        if(bond(iv(i-1))==0) then
          clusterL(k1)=i
        end if
      end if
      if(bond(iv(i-1))==0 .and. bond(i)==1 ) then
        k1=k1+1
        clusterL(k1)=i
      end if
      if(bond(iv(i+1))==0 .and. bond(i)==1 ) then
        k2=k2+1
        clusterR(k2)=i
      end if
    end do
    if(k1.ne.k2) then
      write(*,*) "WARNING!"
    end if
    if(clusterL(1)>clusterR(1) ) then
      do i=1,k2-1
        clusterRaux(i)=clusterR(iv(i+1) )
      end do
      clusterRaux(k2)=clusterR(1)
      clusterR(:)=0._dp
      clusterR(:)=clusterRaux(:)
    end if
    do i=1,k1
      if(clusterL(1).le. clusterR(1) ) then
        write(*,*) clusterL(i),clusterR(i)
        do j=clusterL(i),clusterR(i)
          Sx(iv(j+1) )=Sx2(iv(j+1) )
          Sy(iv(j+1) )=Sy2(iv(j+1) )
        end do
      else
        write(*,*) "Wtf?"
      end if
    end do
    deallocate(bond,Sx2,Sy2,clusterL,clusterR,clusterRaux)
  end subroutine Clustert

   subroutine mean_0(x,y)
    real(dp), dimension(Nmsrs), intent(in) :: x
    real(dp), intent(out) :: y
    integer(i4) :: k
    y=0._dp
    do k=1,Nmsrs
      y=y+x(k)
    end do
    y=y/real(Nmsrs,dp)
  end subroutine mean_0

  subroutine standard_error(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp), intent(out) :: deltay
    real(dp) :: variance
    integer(i4) :: k,N
    N=size(x)
    deltay=0._dp
    variance=0._dp
    do k=1,N
      variance=variance+(x(k) -y)**2
    end do
    variance=variance/real(N-1,dp)
    deltay=Sqrt(variance/real(N,dp))
  end subroutine standard_error

  subroutine jackknife(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp), intent(out) :: deltay
    real(dp) :: jackk
    real(dp),dimension(Mbins) :: xmean
    integer(i4) :: k,N,i
      N=size(x)
      deltay=0._dp
      jackk=0._dp
      xmean=0._dp
      do i=1,Mbins
        do k=1,N
          if(k .le. (i-1)*N/Mbins) then
            xmean(i)=xmean(i)+x(k)
          else if(k > i*N/Mbins) then
            xmean(i)=xmean(i)+x(k)
          end if
        end do
        xmean(i)=xmean(i)/(real(N,dp) -real(N/Mbins,dp))
      end do
      do k=1,Mbins
        jackk=jackk+(xmean(k)-y )**2
      end do
      deltay=Sqrt(real(Mbins-1,dp)*jackk/real(Mbins,dp))
  end subroutine jackknife

  subroutine mean_scalar(x,y,deltay)
    real(dp), dimension(Nmsrs), intent(in) :: x
    real(dp), intent(out) :: y,deltay
    call mean_0(x,y)
    !call standard_error(x,y,deltay)
    call jackknife(x,y,deltay)
  end subroutine mean_scalar


end module statistics
