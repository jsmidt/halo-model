module hm_init
use Precision
use CAMB
implicit none

! Use camb derived type for halo-model parameters HP.
!integer,parameter :: dp = selected_real_kind(p=15,r=307)

! The global halo-model parameter structure. Consistant with CAMB's. 
type, extends(CAMBdata) :: HMdata
   real(dl) :: rho_c
   real(dl) :: r_s
   real(dl) :: z
   real(dl), allocatable,dimension(:) :: k,m,nu_fnu,nu,sig_2,bias_1,bias_2
   real(sp), allocatable,dimension(:) :: Pk
   character*8 :: run_name 
end type HMdata

type(HMdata) :: hm

contains

subroutine init_params()
    !type(CAMBdata) :: hm
    ! Get Default Cosmology Parameters. Match to WMAP 7.
    call CAMB_SetDefParams(hm%Params)

    hm%run_name = 'test'
    hm%z = 0.0
    hm%Params%WantTransfer= .true.
    hm%Params%WantCls= .false.
    hm%Params%Transfer%redshifts=hm%z
    hm%Params%InitPower%ScalarPowerAmp(1) = 2.50e-9
    hm%Params%omegab  = 0.045
    hm%Params%omegac  = 0.227
    hm%Params%omegav  = 0.728
    hm%Params%omegan  = 0
    hm%Params%H0      = 70.2
    hm%Params%YHe     = 0.24
    hm%Params%Num_Nu_massless = 3.04
    hm%Params%Num_Nu_massive  = 0


    ! To get the critical density rho_c we need to calculate 3H^2/(8\piG).  
    ! We want units of M_\sun/Mpc^3. We know H0 has units of km/s/Mpc and 
    ! G = 4.302e-9 Mpc(km/s)^2/M_\sun So that 3H^2/(8\piG) has the right units.
    hm%rho_c = 3.0d0*hm%Params%H0**2/8.0d0/pi/4.302e-9
end subroutine init_params

end module hm_init
