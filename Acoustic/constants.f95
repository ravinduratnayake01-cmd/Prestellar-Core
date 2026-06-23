    MODULE CONSTANTS

    IMPLICIT NONE

    INTEGER, PUBLIC, PARAMETER :: dp = SELECTED_REAL_KIND(P=8)

!   Math constants
    REAL (KIND=dp), PUBLIC, PARAMETER :: pi    = 3.1415926535897932384_dp       ! Pi
    REAL (KIND=dp), PUBLIC, PARAMETER :: sqpi  = 1.77245385090551602_dp         ! Sqrt(pi)
    REAL (KIND=dp), PUBLIC, PARAMETER :: twopi = 6.2831853071795864769_dp       ! 2*Pi

!   Physical constants (CGS)
    REAL (KIND=dp), PUBLIC, PARAMETER :: c_lum    = 2.99791458e+10_dp           ! Speed of light       [cm/s]
    REAL (KIND=dp), PUBLIC, PARAMETER :: h_planck = 6.62606876e-27_dp           ! Planck constant      [erg.s]
    REAL (KIND=dp), PUBLIC, PARAMETER :: kb       = 1.3806503e-16_dp            ! Boltzmann constant   [erg/K]
    REAL (KIND=dp), PUBLIC, PARAMETER :: NAvogadro= 6.022e23_dp                 ! Avogadro number      [mol-1]
    REAL (KIND=dp), PUBLIC, PARAMETER :: mH       = 1.67262178e-24_dp           ! Proton mass          [g]
    REAL (KIND=dp), PUBLIC, PARAMETER :: GNewton  = 6.674e-8_dp                 ! Gravitational const  [dyn cm² g-2]

!   Conversion factors
    REAL (KIND=dp), PUBLIC, PARAMETER :: G0             = 5.6e-14_dp            ! G0
    REAL (KIND=dp), PUBLIC, PARAMETER :: eV_to_erg      = 1.602176487e-12_dp    ! eV -> erg
    REAL (KIND=dp), PUBLIC, PARAMETER :: cm_to_A        = 1.0e+8_dp             ! cm -> Angstrom
    REAL (KIND=dp), PUBLIC, PARAMETER :: cmm1_to_eV     = 1.2398389635e-4_dp    ! cm-1 -> eV
    REAL (KIND=dp), PUBLIC, PARAMETER :: A_to_eV        = h_planck * c_lum * cm_to_A / eV_to_erg ! Angstrom -> eV
    REAL (KIND=dp), PUBLIC, PARAMETER :: parsec_to_cm   = 3.08567758e18_dp      ! [cm/pc]
    REAL (KIND=dp), PUBLIC, PARAMETER :: AU_to_cm       = 1.495978707e13_dp     ! [cm/AU]
    REAL (KIND=dp), PUBLIC, PARAMETER :: yr_to_s        = 3.15576e7_dp          ! [s/yr]

    REAL (KIND=dp), PUBLIC, PARAMETER :: DEGREE_TO_RADIAN = pi / 180.0_dp

!   Astrophysical parameters
    REAL (KIND=dp), PUBLIC, PARAMETER :: muH    = 1.4_dp                        ! Mean particle mass per H nucleus [dimensionless]
    REAL (KIND=dp), PUBLIC, PARAMETER :: fractal= 0.3_dp                        ! Effective fractal optical depth fraction

    END MODULE CONSTANTS
