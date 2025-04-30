module functions
    use iso_fortran_env, only : dp => real64, i4 => int32
    use parameters
    implicit none

contains

  function iv(i)
    integer(i4), intent(in) :: i
    integer(i4) :: iv
    if(i==L+1) then
      iv=1
    else if(i==0) then
      iv=L
    else
      iv=i
    end if
  end function iv

  function Hamilt(Sx,Sy)
    real(dp), dimension(:), intent(in) :: Sx,Sy
    real(dp) :: Hamilt
    integer(i4) :: i,Narr
    Narr=size(Sx)
    Hamilt=0._dp
    do i=1,Narr
      Hamilt=Hamilt-Sx(i)*Sx(iv(i+1) )-Sy(i)*Sy(iv(i+1) )
    end do
    Hamilt=Hamilt
  end function Hamilt

  function DeltaE(Sx,Sy,i,Sx2,Sy2)
    real(dp), dimension(:), intent(in) :: Sx,Sy
    integer(i4), intent(in) :: i
    real(dp), intent(in) :: Sx2,Sy2
    real(dp) :: DeltaEx,DeltaEy,DeltaE
      DeltaEx=(Sx(iv(i+1) )+Sx(iv(i-1) ) )*(Sx(i)-Sx2 )
      DeltaEy=(Sy(iv(i+1) )+Sy(iv(i-1) ) )*(Sy(i)-Sy2 )
      DeltaE=(DeltaEx+DeltaEy)
  end function DeltaE

  function Deltah(Sx,Sy,i,Sx2,Sy2)
    real(dp), dimension(:), intent(in) :: Sx,Sy
    integer(i4), intent(in) :: i
    real(dp), intent(in) :: Sx2,Sy2
    real(dp) :: Deltah
    Deltah=(Sx(i)*(Sx(iv(i+1))-Sx2 )+Sy(i)*(Sy(iv(i+1))-Sy2 ) )
  end function Deltah

  function Magnet2(Sx,Sy)
    real(dp), dimension(:), intent(in) :: Sx,Sy
    real(dp) :: Magnet2,a,b
    integer(i4) :: i,Narr
    Narr=size(Sx)
    Magnet2=0._dp
    a=0._dp
    b=0._dp
    do i=1,Narr
      a=a+Sx(i)
      b=b+Sy(i)
    end do
    Magnet2=(a**2+b**2)/real(Narr,dp)
  end function Magnet2

  function top_charge(Sx,Sy)
    real(dp), dimension(L), intent(in) :: Sx,Sy
    real(dp), dimension(L) :: phi(:)
    real(dp) :: dphi
    integer(i4) :: i
    do i=1,L
      phi(i)=atan2(Sy(i),Sx(i))
    end do
    do i=1,L
      dphi=dphi+phi(iv(i+1))-phi(i)
    end do
  end function top_charge


end module functions
