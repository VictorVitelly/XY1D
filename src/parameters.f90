module parameters

    use iso_fortran_env, only : dp => real64, i4 => int32

    integer(i4), parameter :: L=100
    real(dp), parameter :: Pi = 4._dp*atan(1.0_dp)
    integer(i4), parameter :: thermalization=500,eachsweep=50,Nmsrs=100000,bins=11,Mbins=10
    integer(i4) :: sweeps=thermalization+eachsweep*Nmsrs
    real(dp) :: dtheta=0.1_dp*Pi

    real(dp), parameter :: maxx=5.5_dp, minn=-5.5_dp
    real(dp) :: binwidth=(maxx-minn)/real(bins,dp)

    real :: starting,ending

end module parameters
