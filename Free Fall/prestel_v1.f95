MODULE PRESTEL

   USE CONSTANTS
   IMPLICIT NONE

   REAL(KIND=dp), PARAMETER :: DENSITY_FLOOR = 1.0e-40_dp   ! minimum density

CONTAINS

!==============================================================================
   SUBROUTINE PRESTEL_STEP(rho_in, v_in, r, cs, dt, Npts, rho_out, v_out)

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: Npts
      REAL(KIND=dp),  INTENT(IN)  :: rho_in(Npts), v_in(Npts), r(Npts)
      REAL(KIND=dp),  INTENT(IN)  :: cs, dt
      REAL(KIND=dp),  INTENT(OUT) :: rho_out(Npts), v_out(Npts)

      ! Local work arrays
      REAL(KIND=dp) :: flux(Npts)         ! r^2 * rho * v
      REAL(KIND=dp) :: dflux_dr(Npts)     ! d(r^2 rho v)/dr
      REAL(KIND=dp) :: drho_dr(Npts)
      REAL(KIND=dp) :: dv_dr(Npts)
      REAL(KIND=dp) :: grav(Npts)         ! gravitational acceleration
      REAL(KIND=dp) :: drho_dt(Npts)
      REAL(KIND=dp) :: dv_dt(Npts)
      INTEGER :: i

      ! Build spherical flux  f = r^2 * rho * v 
      DO i = 1, Npts
         flux(i) = r(i)**2 * rho_in(i) * v_in(i)
      END DO

      ! Centred gradients 
      CALL GRADIENT_SUB(flux, r, Npts, dflux_dr)
      CALL GRADIENT_SUB(rho_in, r, Npts, drho_dr)
      CALL GRADIENT_SUB(v_in, r, Npts, dv_dr)

      ! Gravitational acceleration g(r) = G*M(<r)/r^2 
      CALL GRAV_FORCE_SUB(rho_in, r, Npts, grav)

      ! Time derivatives (explicit Euler RHS) 
      DO i = 1, Npts
         ! Continuity: d(rho)/dt = -(1/r^2) * d(r^2*rho*v)/dr
         drho_dt(i) = -(1.0_dp / r(i)**2) * dflux_dr(i)

         ! Momentum: d(v)/dt = -v*dv/dr - (cs^2/rho)*drho/dr - g
         dv_dt(i) = -v_in(i) * dv_dr(i) &
                    - (cs**2 / rho_in(i)) * drho_dr(i) &
                    - grav(i)
      END DO

      ! Explicit Euler step 
      rho_out = rho_in + drho_dt * dt
      v_out = v_in + dv_dt   * dt

      ! Reflective inner boundary (no inflow at r=0) 
      rho_out(1) = rho_out(2)
      v_out(1)   = -v_out(2)          ! antisymmetric velocity at origin

      ! Outflow / zero-gradient outer boundary 
      rho_out(Npts) = rho_out(Npts-1)
      v_out(Npts) = v_out(Npts-1)

      ! Density floor 
      WHERE (rho_out <= 0.0_dp) rho_out = DENSITY_FLOOR

   END SUBROUTINE PRESTEL_STEP

!==============================================================================
   SUBROUTINE GRADIENT_SUB(F, r, Npts, dFdr)

      IMPLICIT NONE

      INTEGER, INTENT(IN)        :: Npts
      REAL(KIND=dp), INTENT(IN)  :: F(Npts), r(Npts)
      REAL(KIND=dp), INTENT(OUT) :: dFdr(Npts)

      INTEGER :: i

      ! Forward difference at left boundary
      dFdr(1) = (F(2) - F(1)) / (r(2) - r(1))

      ! Centred differences for interior points
      DO i = 2, Npts - 1
         dFdr(i) = (F(i+1) - F(i-1)) / (r(i+1) - r(i-1))
      END DO

      ! Backward difference at right boundary
      dFdr(Npts) = (F(Npts) - F(Npts-1)) / (r(Npts) - r(Npts-1))

   END SUBROUTINE GRADIENT_SUB

!==============================================================================
   SUBROUTINE GRAV_FORCE_SUB(rho, r, Npts, grav)

      IMPLICIT NONE

      INTEGER,        INTENT(IN)  :: Npts
      REAL(KIND=dp),  INTENT(IN)  :: rho(Npts), r(Npts)
      REAL(KIND=dp),  INTENT(OUT) :: grav(Npts)

      INTEGER :: i
      REAL(KIND=dp) :: mass_enc, dr

      mass_enc = 0.0_dp
      grav(1) = 0.0_dp          ! No mass enclosed at innermost cell

      DO i = 2, Npts
         dr = r(i) - r(i-1)
         mass_enc = mass_enc + 0.5_dp * 4.0_dp * pi * &
                    (rho(i) * r(i)**2 + rho(i-1) * r(i-1)**2) * dr
         grav(i)  = GNewton * mass_enc / r(i)**2
      END DO

   END SUBROUTINE GRAV_FORCE_SUB


!==============================================================================
   FUNCTION TOTAL_MASS(rho, r, Npts) RESULT(mass)

      IMPLICIT NONE

      INTEGER,        INTENT(IN) :: Npts
      REAL(KIND=dp),  INTENT(IN) :: rho(Npts), r(Npts)
      REAL(KIND=dp) :: mass

      INTEGER :: i

      mass = 0.0_dp
      DO i = 2, Npts
         mass = mass + 0.5_dp * 4.0_dp * pi * &
                (rho(i) * r(i)**2 + rho(i-1) * r(i-1)**2) * (r(i) - r(i-1))
      END DO

   END FUNCTION TOTAL_MASS


!==============================================================================
   SUBROUTINE RHO_TO_NH(rho, Npts, nH)

      IMPLICIT NONE

      INTEGER,        INTENT(IN)  :: Npts
      REAL(KIND=dp),  INTENT(IN)  :: rho(Npts)
      REAL(KIND=dp),  INTENT(OUT) :: nH(Npts)

      nH = rho / (muH * mH)

   END SUBROUTINE RHO_TO_NH


!  Convert hydrogen number density nH [cm^-3] to mass density rho [g/cm^3]
!==============================================================================
   SUBROUTINE NH_TO_RHO(nH, Npts, rho)

      IMPLICIT NONE

      INTEGER,        INTENT(IN)  :: Npts
      REAL(KIND=dp),  INTENT(IN)  :: nH(Npts)
      REAL(KIND=dp),  INTENT(OUT) :: rho(Npts)

      rho = nH * muH * mH

   END SUBROUTINE NH_TO_RHO


!  Mass accretion rate at grid index idx
!==============================================================================
   FUNCTION ACCRETION_RATE(rho, v, r, Npts, idx) RESULT(mdot)

      IMPLICIT NONE

      INTEGER,        INTENT(IN) :: Npts, idx
      REAL(KIND=dp),  INTENT(IN) :: rho(Npts), v(Npts), r(Npts)
      REAL(KIND=dp) :: mdot

      IF (idx < 1 .OR. idx > Npts) THEN
         WRITE(*,'(A,I6,A,I6,A)') &
            'ERROR ACCRETION_RATE: index ', idx, ' out of range [1,', Npts, ']'
         STOP
      END IF

      mdot = 4.0_dp * pi * r(idx)**2 * rho(idx) * v(idx)

   END FUNCTION ACCRETION_RATE

!  CFL_CHECK
!==============================================================================
   SUBROUTINE CFL_CHECK(v, r, cs, Npts, CFL_factor, dt_user, dt_cfl, ok)

      IMPLICIT NONE

      INTEGER,        INTENT(IN)  :: Npts
      REAL(KIND=dp),  INTENT(IN)  :: v(Npts), r(Npts)
      REAL(KIND=dp),  INTENT(IN)  :: cs, CFL_factor, dt_user
      REAL(KIND=dp),  INTENT(OUT) :: dt_cfl
      LOGICAL,        INTENT(OUT) :: ok

      INTEGER :: i
      REAL(KIND=dp) :: dr, sig, min_ratio

      min_ratio = HUGE(1.0_dp)
      DO i = 1, Npts - 1
         dr = r(i+1) - r(i)
         sig = ABS(v(i)) + cs
         IF (sig > 0.0_dp) min_ratio = MIN(min_ratio, dr / sig)
      END DO

      dt_cfl = CFL_factor * min_ratio
      ok = (dt_user <= dt_cfl)

   END SUBROUTINE CFL_CHECK


!==============================================================================
   FUNCTION SOUND_SPEED(T) RESULT(cs)

      IMPLICIT NONE

      REAL(KIND=dp), INTENT(IN) :: T
      REAL(KIND=dp) :: cs

      cs = SQRT(kb * T / (muH * mH))

   END FUNCTION SOUND_SPEED


!  Build a 1D radial grid.
!==============================================================================
   SUBROUTINE RADIAL_GRID(Npts, Rmin, Rmax, gridtype, r)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: Npts
      REAL(KIND=dp), INTENT(IN) :: Rmin, Rmax
      CHARACTER(LEN=3), INTENT(IN) :: gridtype
      REAL(KIND=dp), INTENT(OUT) :: r(Npts)

      INTEGER :: i
      REAL(KIND=dp) :: c

      SELECT CASE (TRIM(gridtype))

         CASE ('lin', 'LIN')
            c = (Rmax - Rmin) / REAL(Npts - 1, dp)
            DO i = 1, Npts
               r(i) = Rmin + c * REAL(i - 1, dp)
            END DO

         CASE ('log', 'LOG')
            c = LOG10(Rmax / Rmin) / REAL(Npts - 1, dp)
            DO i = 1, Npts
               r(i) = Rmin * 10.0_dp**(c * REAL(i - 1, dp))
            END DO

         CASE DEFAULT
            WRITE(*,'(A,A)') 'ERROR RADIAL_GRID: unknown grid type: ', gridtype
            STOP

      END SELECT

   END SUBROUTINE RADIAL_GRID

END MODULE PRESTEL

