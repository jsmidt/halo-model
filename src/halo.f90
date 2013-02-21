!
!  halo.f90 - calculates the n-halo terms for the halo model.
!
!    Halo model calculations of n-point functions involve integrals 
!    over weighted mass functions and Fourier transforms of the
!    underlying halo distributions.  The contributions of these n-halo
!    terms are performed in this module. For more informations of each
!    of these n-halo terms for the power, bispectrum and trispectrum 
!    see eqs. 95-98 of Cooray & Sheth (2002) astro-ph/0206508.
!
!  Version history:
!    v0.6  Joseph Smidt 05/31/2012 - Initial 1 and 2 halo terms for P(k)
!
!  See also: massfunc.f90 and hm_init.f90
!
!  Copyright (c) 2012 Joseph Smidt. See accompanying BSD license.
!



module halo
use massfunc
implicit none


contains

! Get 1-Halo Term P^{1h}(k).  See eq. 95 of Cooray and Sheth (2002).
real(dl) function P1h(kg)
    real(dl), intent(in) :: kg
    real(dl) :: rombint,tol2

    ! Calculate the integral over integrand P1hi over \nu. 
    ! M_{0,2} of eq. 95 of Cooray and Sheth (2002).
    CALL qromb(P1hi,dlog(minval(hm%nu_m)),dlog(maxval(hm%nu_m)),P1h,kg)

end function P1h
real(dl) function P1hi(lnnu,k)
    real(dl),intent(in) :: lnnu
    real(dl) :: m,k,inu_fnu,omegam,omegal,omegamz,nu,iukm

    ! Get \Omega_m(z) to get mass densty \rho_m = \rho_c\Omega_m(z)
    omegam = hm%Params%omegab + hm%Params%omegac
    omegal = hm%Params%omegav
    omegamz = omegam*(1.0+hm%z)**3/(omegam*(1.0+hm%z)**3+omegal)
    nu = exp(lnnu)

    ! Interpolate to get values at each \nu.
    m = interpf(dlog(hm%nu_m),hm%m,lnnu)
    call splint(dlog(hm%nu_m),hm%nu_fnum,hm%sp_nu_fnum,size(hm%nu_m), &
            lnnu,inu_fnu)
    call splin2(dlog(hm%k),dlog(hm%m),hm%ukm,hm%sp_ukm,hm%mpts,hm%mpts, &
            dlog(k),dlog(m),iukm)

    ! Integrand is \nuf_{\nu}u(k|m)^2
    P1hi = m/hm%rho_c/omegamz*inu_fnu*iukm**2

end function P1hi

! Get 2-Halo Term P^{2h}(k).  See eq. 95 of Cooray and Sheth (2002).
real(dl) function P2h(kg)
    real(dl), intent(in) :: kg

    ! Calculate the integral over integrand P2hi over \nu. 
    ! M_{1,1} of eq. 95 of Cooray and Sheth (2002).
    CALL qromb(P2hi,log(minval(hm%nu_m)),dlog(maxval(hm%nu_m)),P2h,kg)

    ! Multiply M_{1,1}^2 by P(k)^{lin}
    P2h = P2h**2*interpf(dlog(hm%k),dble(hm%Pk),log(kg))

end function P2h
real(dl) function P2hi(lnnu,k)
    real(dl), intent(in) :: lnnu,k
    real(dl) :: m,inu_fnu,ibias_1,nu,iukm
    nu = exp(lnnu)

    ! Interpolate to get values at each \nu.
    m = interpf(dlog(hm%nu_m),hm%m,lnnu)
    call splint(dlog(hm%nu_m),hm%nu_fnum,hm%sp_nu_fnum,size(hm%nu_m),lnnu, &
            inu_fnu)
    call splint(dlog(hm%nu_m),hm%bias_1,hm%sp_bias_1,size(hm%nu_m),lnnu,&
            ibias_1)
    call splin2(dlog(hm%k),dlog(hm%m),hm%ukm,hm%sp_ukm,hm%mpts,hm%mpts, &
            dlog(k),dlog(m),iukm)

    ! Integrand is \nuf_{\nu}b_1{k)u(k|m)
    P2hi = ibias_1*inu_fnu*iukm

end function P2hi

u
! End module.
end module halo
