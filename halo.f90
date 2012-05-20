module halo
use CAMB
implicit none

integer, parameter::dp = kind(1.0d0)
real(dp), parameter:: kmin = 2.0e-4
real(dp), parameter:: kmax = 9.610024d0
real(dp), parameter:: mmin = 1e8
real(dp), parameter:: mmax = 1e13
real(dp), parameter:: dlnk = 0.07d0
integer, parameter:: mpts = 155
real(dp), parameter:: rho_bar = 1.0d0
real(dp) :: R,z

contains

! Tophat window function. See eq. 58 of astro-ph/0206508.
elemental real(dp) function w_top(x)
    real(dp),intent(in) :: x
    w_top = (3.0d0/x**3)*(sin(x)-x*cos(x))
end function w_top


! Variance in the initial density fluctuation (\sigma^2(m))
! Equation 58 of astro-ph/0206508.
real(dp) function sig_2(m)
    real(dp),intent(in) :: m
    real(dp) :: a,b,rombint
    a = kmin
    b = kmax
    R = (3.0d0*m/4.0d0/pi/rho_bar)**0.3333333333d0
    sig_2 = rombint(sig_2int,a,b,tol)
end function sig_2
real(dp) function sig_2int(k)
    real(dp),intent(in) :: k
    sig_2int = (k**2/2.0d0/pi**2)*P_lin(k)*w_top(k*R)**2
end function sig_2int
real(dp) function P_lin(k)
real(dp), intent(in) :: k
P_lin = 2
end function P_lin

! Calculate nu for mass m and redshift z.  Equation 57 of astro-ph/0206508.
real(dp) function nu(m)
    real(dp), intent(in) :: m
    nu = 1.686470199841145d0*(1.0d0+z)*sig_2(m)
end function nu

! Calculate f_nu for mass m and redshift z. Equation 59 of astro-ph/0206508.
real(dp) function fnu(m)
    real(dp), intent(in) :: m
    real(dp):: A,q
    A = 0.3222d0
    q = 0.75d0
    fnu=A*(1.0d0+(q*nu(m))**0.3)*(q*nu(m)/2.0d0/pi)**0.5*exp(-q*nu(m)/2.0d0)/nu(m)
end function fnu

! Calculate the Fourier transform of the dark matter distribution u(k|m)
! Equation 81 & 82 of astro-ph/0206508
real(dp) function ukm(k,m)
    real(dp), intent(in) :: k
    real(dp) :: m
    real(dp) :: ps, rs,c 
    ps = 1.0
    rs = 1.0
    c = 1.0
    ukm = 4.0d0*pi*ps*rs**3/m*(sin(k*rs)*(Si((1.0d0+C)*k*rs)-Si(k*rs)) &
        - sin(c*k*rs)/((1.0d0+c)*k*rs)+cos(k*rs)*(Ci((1.0d0+C)*k*rs)-Ci(k*rs)))
end function ukm
real(dp) function Ci(x)
   real(dp), intent(in) :: x
   real(dp) :: rombint
   Ci = rombint(Cii,x,kmax,tol)
end function Ci
real(dp) function Si(x)
   real(dp), intent(in) :: x
   real(dp) :: rombint
   Si = rombint(Sii,x,kmax,tol)
end function Si
real(dp) function Cii(t)
   real(dp), intent(in) :: t
   Cii = -cos(t)/t
end function Cii
real(dp) function Sii(t)
   real(dp), intent(in) :: t
   Sii = sin(t)/t
end function Sii


! Get linear power spectrum Pk at k from CAMB. 
subroutine linear_pk(k,Pk)
   real, allocatable,dimension(:),intent(inout) :: k,Pk
   integer :: i,error
   type(CAMBdata) :: P

   !allocate(k(mpts),Pk(mpts))
   do i=1,mpts
     k(i)=exp(log(kmin)+dlnk*(i-1))
   end do

   ! Get Default Cosmology Parameters.
   call CAMB_SetDefParams(P%Params)
   P%Params%WantTransfer= .true.
   P%Params%Transfer%redshifts=z

   ! Get transfer functions and power spectra Pk
   call CAMB_GetTransfers(P%Params, P, error)
   call Transfer_GetMatterPower(P%MTrans,Pk,1,1,real(kmin),real(dlnk),mpts)
   Pk = Pk/3.0e8

end subroutine linear_pk

! Get 1-Halo Term
real(dp) function P1h(k)
    real(dp), intent(in) :: k
    real(dp) :: rombint,tol2
    tol2 = 1e-2
    P1h = rombint(P1hi,mmin,mmax,tol2)
end function P1h
real(dp) function P1hi(m)
    real(dp) :: m
    real(dp) :: k
    k = 1.0
    P1hi = fnu(m)*(m/rho_bar)*ukm(k,m)**2
end function P1hi

end module halo
