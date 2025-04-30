module parameters

    use iso_fortran_env, only : dp => real64, i4 => int32

    integer(i4), parameter :: L=50
    real(dp), parameter :: Pi = 4._dp*atan(1.0_dp)
    integer(i4), parameter :: thermalization=10000,eachsweep=50,Nmsrs=10000,bins=101,Mbins=10
    integer(i4) :: sweeps=thermalization+eachsweep*Nmsrs
    real(dp) :: dtheta=0.2_dp*Pi

    !real(dp), parameter :: maxx=2.0_dp, minn=-2.0_dp
    !real(dp) :: binwidth=(maxx-minn)/real(bins,dp)

    real :: starting,ending

end module parameters
