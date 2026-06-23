!===============================================================================

!  Jeans instability analysis.


!    omega^2 = cs^2 * k^2 - 4*pi*G*rho0

!===============================================================================

   MODULE JEANS_ANALYSIS

   USE CONSTANTS

   IMPLICIT NONE

   CONTAINS

!===============================================================================
!  Critical wavenumber kJ where omega=0

!  kJ = sqrt(4*pi*G*rho0) / cs

!===============================================================================
   FUNCTION JEANS_WAVENUMBER(rho0, cs)

   IMPLICIT NONE
   REAL(KIND=dp), INTENT(IN) :: rho0, cs
   REAL(KIND=dp) :: JEANS_WAVENUMBER

   JEANS_WAVENUMBER = SQRT(4.0_dp * pi * GNewton * rho0) / cs

   END FUNCTION JEANS_WAVENUMBER
!===============================================================================


!===============================================================================
! Critical wavelength lambdaJ

!  lambdaJ = 2*pi / kJ = cs * sqrt(pi / (G*rho0))
!===============================================================================
   FUNCTION JEANS_LENGTH(rho0, cs)

   IMPLICIT NONE
   REAL(KIND=dp), INTENT(IN) :: rho0, cs
   REAL(KIND=dp) :: JEANS_LENGTH

   JEANS_LENGTH = cs * SQRT(pi / (GNewton * rho0))

   END FUNCTION JEANS_LENGTH
!===============================================================================


!===============================================================================
!  Jeans mass MJ = (4/3)*pi*(lambdaJ/2)^3 * rho0

!  The mass contained in a sphere of diameter lambdaJ.

!===============================================================================
   FUNCTION JEANS_MASS(rho0, cs)

   IMPLICIT NONE
   REAL(KIND=dp), INTENT(IN) :: rho0, cs
   REAL(KIND=dp) :: JEANS_MASS
   REAL(KIND=dp) :: lJ

   lJ = JEANS_LENGTH(rho0, cs)
   JEANS_MASS = (4.0_dp / 3.0_dp) * pi * (lJ / 2.0_dp)**3 * rho0

   END FUNCTION JEANS_MASS
!===============================================================================


!===============================================================================
!  omega^2 from the dispersion relation for given k

!  omega^2 = cs^2 * k^2 - 4*pi*G*rho0

!===============================================================================
   FUNCTION DISPERSION_OMEGA2(k, rho0, cs)

   IMPLICIT NONE
   REAL(KIND=dp), INTENT(IN) :: k, rho0, cs
   REAL(KIND=dp) :: DISPERSION_OMEGA2

   DISPERSION_OMEGA2 = cs**2 * k**2 - 4.0_dp * pi * GNewton * rho0

   END FUNCTION DISPERSION_OMEGA2
!===============================================================================


!===============================================================================
!  Exponential growth rate sigma for unstable perturbations

!  For omega^2 < 0, the perturbation grows as exp(sigma * t) where:
!    sigma = sqrt(4*pi*G*rho0 - cs^2*k^2)

!===============================================================================
   FUNCTION GROWTH_RATE(k, rho0, cs)

   IMPLICIT NONE
   REAL(KIND=dp), INTENT(IN) :: k, rho0, cs
   REAL(KIND=dp) :: GROWTH_RATE
   REAL(KIND=dp) :: omega2

   omega2 = DISPERSION_OMEGA2(k, rho0, cs)

   IF (omega2 < 0.0_dp) THEN
      GROWTH_RATE = SQRT(-omega2)
   ELSE
      GROWTH_RATE = 0.0_dp
   END IF

   END FUNCTION GROWTH_RATE
!===============================================================================


!===============================================================================
!  Characteristic collapse timescale ~ 1/sigma )

!  The maximum growth rate is at k=0:  sigma_max = sqrt(4*pi*G*rho0), free-fall timescale

!  t_J = 1 / sqrt(4*pi*G*rho0)
!===============================================================================
   FUNCTION JEANS_TIME(rho0)

   IMPLICIT NONE
   REAL(KIND=dp), INTENT(IN) :: rho0
   REAL(KIND=dp) :: JEANS_TIME

   JEANS_TIME = 1.0_dp / SQRT(4.0_dp * pi * GNewton * rho0)

   END FUNCTION JEANS_TIME
!===============================================================================



   END MODULE JEANS_ANALYSIS
