!===========================
   MODULE PRESTEL
!===========================

!===============================================================================

   USE CONSTANTS

   IMPLICIT NONE

!===============================================================================

   CONTAINS

!===============================================================================
!  Build radial grid, linear and logarithmic

!    Npts - number of grid points
!    Rmin - inner boundary [cm]
!    Rmax - outer boundary [cm]
!    gridtype - 'lin' for linear, 'log' for logarithmic

!===============================================================================
   FUNCTION RADIAL_GRID(Npts, Rmin, Rmax, gridtype)

   IMPLICIT NONE

   INTEGER,          INTENT(IN) :: Npts
   REAL (KIND=dp),   INTENT(IN) :: Rmin, Rmax
   CHARACTER(LEN=3), INTENT(IN) :: gridtype
   REAL (KIND=dp), DIMENSION(Npts) :: RADIAL_GRID

   INTEGER :: i
   REAL (KIND=dp) :: c

   SELECT CASE (TRIM(gridtype))

      CASE ('lin', 'LIN')
         c = (Rmax - Rmin) / REAL(Npts - 1, dp)
         DO i = 1, Npts
            RADIAL_GRID(i) = Rmin + c * REAL(i - 1, dp)
         END DO

      CASE ('log', 'LOG')
         ! Logarithmically spacing
         ! Rmin > 0
         IF (Rmin <= 0.0_dp) THEN
            WRITE(*,*) 'Rmin must be > 0 for log grid.'
            STOP
         END IF
         c = LOG10(Rmax / Rmin) / REAL(Npts - 1, dp)
         DO i = 1, Npts
            RADIAL_GRID(i) = Rmin * 10.0_dp ** (c * REAL(i - 1, dp))
         END DO

      CASE DEFAULT
         WRITE(*,*) 'Unknown grid type "', TRIM(gridtype), '"'
         WRITE(*,*) '  Supported types: "lin", "log"'
         STOP

   END SELECT

   END FUNCTION RADIAL_GRID
!===============================================================================


!===============================================================================
!  Isothermal sound speed
!  cs = sqrt( kB * T / (muH * mH) )

!    T  - temperature [K]
!===============================================================================
   FUNCTION SOUND_SPEED(T)

   IMPLICIT NONE

   REAL (KIND=dp), INTENT(IN) :: T
   REAL (KIND=dp) :: SOUND_SPEED

   SOUND_SPEED = SQRT(kb * T / (muH * mH))

   END FUNCTION SOUND_SPEED
!===============================================================================


!===============================================================================
!  Centred finite difference gradient of a 1D field F on grid r

!  Uses centred differences for interior points:
!    dF/dr(i) = [ F(i+1) - F(i-1) ] / [ r(i+1) - r(i-1) ]
!  One-sided differences at the boundaries:
!    dF/dr(1) = [ F(2)    - F(1)    ] / [ r(2)    - r(1)    ]   (forward)
!    dF/dr(Npts) = [ F(Npts) - F(Npts-1) ] / [ r(Npts) - r(Npts-1) ] (backward)

!    F(Npts) - field to differentiate
!    r(Npts) - radial grid [cm]
!    Npts  - number of grid points

!===============================================================================
   FUNCTION GRADIENT(F, r, Npts)

   IMPLICIT NONE

   INTEGER,        INTENT(IN) :: Npts
   REAL (KIND=dp), INTENT(IN) :: F(Npts), r(Npts)
   REAL (KIND=dp) :: GRADIENT(Npts)

   INTEGER :: i

   ! Forward difference at left boundary
   GRADIENT(1) = (F(2) - F(1)) / (r(2) - r(1))

   ! Centred differences for interior
   DO i = 2, Npts - 1
      GRADIENT(i) = (F(i+1) - F(i-1)) / (r(i+1) - r(i-1))
   END DO

   ! Backward difference at right boundary
   GRADIENT(Npts) = (F(Npts) - F(Npts-1)) / (r(Npts) - r(Npts-1))

   END FUNCTION GRADIENT
!===============================================================================


!===============================================================================
!  Compute gravitational acceleration dPhi/dr = G * M(<r) / r^2

!    dPhi/dr(r) = G * M(<r) / r^2
!   M(<r) = 4*pi * integral_0^r rho(r') r'^2 dr'

!    rho(Npts) - mass density profile [g/cm^3]
!    r(Npts)   - radial grid [cm]
!    Npts      - number of grid points
!
!===============================================================================
   FUNCTION GRAV_FORCE(rho, r, Npts)

   IMPLICIT NONE

   INTEGER,        INTENT(IN) :: Npts
   REAL (KIND=dp), INTENT(IN) :: rho(Npts), r(Npts)
   REAL (KIND=dp) :: GRAV_FORCE(Npts)

   INTEGER :: i
   REAL (KIND=dp) :: mass_enclosed, dr

   ! M(<r) by cumulative trapezoidal integration of 4*pi*r^2*rho
   mass_enclosed = 0.0_dp
   GRAV_FORCE(1) = 0.0_dp    ! Since no enclosed mass at innermost point

   DO i = 2, Npts
      dr = r(i) - r(i-1)
      ! Trapezoidal rule, 0.5*(f(i)+f(i-1))*dr,  f = 4*pi*r^2*rho
      mass_enclosed = mass_enclosed + 0.5_dp * 4.0_dp * pi * ( rho(i)   * r(i)**2   + rho(i-1) * r(i-1)**2 ) * dr
      GRAV_FORCE(i) = GNewton * mass_enclosed / r(i)**2
   END DO

   END FUNCTION GRAV_FORCE
!===============================================================================


!===============================================================================
!  Total mass of the core by spherical integration

!  M = 4*pi * integral_0^Rmax rho(r) r^2 dr   trapezoidal

!    rho(Npts) - mass density profile [g/cm^3]
!    r(Npts)   - radial grid [cm]
!    Npts      - number of grid points

!===============================================================================
   FUNCTION TOTAL_MASS(rho, r, Npts)

   IMPLICIT NONE

   INTEGER,        INTENT(IN) :: Npts
   REAL (KIND=dp), INTENT(IN) :: rho(Npts), r(Npts)
   REAL (KIND=dp) :: TOTAL_MASS

   INTEGER :: i

   TOTAL_MASS = 0.0_dp
   DO i = 2, Npts
      TOTAL_MASS = TOTAL_MASS + 0.5_dp * 4.0_dp * pi * ( rho(i)   * r(i)**2   + rho(i-1) * r(i-1)**2 ) * (r(i) - r(i-1))
   END DO

   END FUNCTION TOTAL_MASS
!===============================================================================


!===============================================================================
!  Convert mass density to hydrogen number density

!  nH = rho / (muH * mH)

!  INPUT:
!    rho(Npts) - mass density [g/cm^3]
!    Npts      - number of grid points

!===============================================================================
   FUNCTION RHO_TO_NH(rho, Npts)

   IMPLICIT NONE

   INTEGER,        INTENT(IN) :: Npts
   REAL (KIND=dp), INTENT(IN) :: rho(Npts)
   REAL (KIND=dp) :: RHO_TO_NH(Npts)

   RHO_TO_NH = rho / (muH * mH)

   END FUNCTION RHO_TO_NH
!===============================================================================


!===============================================================================
!  Convert hydrogen number density to mass density

!  rho = nH * muH * mH

!    nH(Npts)  - number density [cm^-3]
!    Npts - number of grid points

!===============================================================================
   FUNCTION NH_TO_RHO(nH, Npts)

   IMPLICIT NONE

   INTEGER,        INTENT(IN) :: Npts
   REAL (KIND=dp), INTENT(IN) :: nH(Npts)
   REAL (KIND=dp) :: NH_TO_RHO(Npts)

   NH_TO_RHO = nH * muH * mH

   END FUNCTION NH_TO_RHO
!===============================================================================


!===============================================================================
!  Mass accretion rate Mdot at a given radial index

!  Mdot = 4 * pi * r^2 * rho * v

!  Positive v mass moving outward
!  For infall, v < 0 and Mdot < 0 mass flowing inward

!    rho(Npts) - mass density profile [g/cm^3]
!    v(Npts) - velocity profile [cm/s]
!    r(Npts)- radial grid [cm]
!    idx - index of the radius at which to evaluate Mdot

!===============================================================================
   FUNCTION ACCRETION_RATE(rho, v, r, Npts, idx)

   IMPLICIT NONE

   INTEGER,        INTENT(IN) :: Npts, idx
   REAL (KIND=dp), INTENT(IN) :: rho(Npts), v(Npts), r(Npts)
   REAL (KIND=dp) :: ACCRETION_RATE

   IF (idx < 1 .OR. idx > Npts) THEN
      WRITE(*,*) 'ERROR in ACCRETION_RATE: index', idx, 'out of range [1,', Npts, ']'
      STOP
   END IF

   ACCRETION_RATE = 4.0_dp * pi * r(idx)**2 * rho(idx) * v(idx)

   END FUNCTION ACCRETION_RATE
!===============================================================================


!===============================================================================
!  Courant-Friedrichs-Lewy timestep condition
!
!  The CFL condition requires that information does not travel more than one grid cell per timestep. The maximum,

!    dt_CFL = CFL_factor * min_i( dr(i) / (|v(i)| + cs) )

!  where CFL_factor < 1 ===> ~ 0.5 is a safety factor

!    v(Npts) - velocity profile [cm/s]
!    r(Npts)- radial grid [cm]
!    cs  - sound speed [cm/s]
!    Npts - number of grid points
!    CFL_factor - safety factor (0 < CFL_factor <= 1, recommend 0.5)

!===============================================================================
   FUNCTION CFL_DT(v, r, cs, Npts, CFL_factor)

   IMPLICIT NONE

   INTEGER,        INTENT(IN) :: Npts
   REAL (KIND=dp), INTENT(IN) :: v(Npts), r(Npts), cs, CFL_factor
   REAL (KIND=dp) :: CFL_DT

   INTEGER :: i
   REAL (KIND=dp) :: dr, signal_speed, min_ratio

   min_ratio = HUGE(1.0_dp)

   DO i = 1, Npts - 1
      dr = r(i+1) - r(i)
      signal_speed = ABS(v(i)) + cs
      IF (signal_speed > 0.0_dp) THEN
         min_ratio = MIN(min_ratio, dr / signal_speed)
      END IF
   END DO

   CFL_DT = CFL_factor * min_ratio

   END FUNCTION CFL_DT
!===============================================================================


!===============================================================================
!  Advance the system by one timestep: Integrate the continuity and Euler equations in 1D spherical symmetry

!    d(rho)/dt = -(1/r^2) * d(r^2 * rho * v)/dr
!    d(v)/dt   = -v * dv/dr - (cs^2/rho) * d(rho)/dr - dPhi/dr
!
!  Centred finite differences and explicit Euler timestepping

!    rho_in(Npts)- density at time t [g/cm^3]
!    v_in(Npts) - velocity at time t [cm/s]
!    r(Npts) - radial grid [cm]  (fixed, does not evolve)
!    cs - isothermal sound speed [cm/s]
!    dt - timestep [s]
!    Npts- number of grid points

!===============================================================================
   SUBROUTINE PRESTEL_STEP_SPHERICAL(rho_in, v_in, r, cs, dt, Npts, use_gravity, rho_out, v_out)

   IMPLICIT NONE

   INTEGER,        INTENT(IN)  :: Npts
   REAL (KIND=dp), INTENT(IN)  :: rho_in(Npts), v_in(Npts), r(Npts)
   REAL (KIND=dp), INTENT(IN)  :: cs, dt
   LOGICAL,        INTENT(IN)  :: use_gravity
   REAL (KIND=dp), INTENT(OUT) :: rho_out(Npts), v_out(Npts)

   ! Local arrays
   REAL (KIND=dp) :: flux(Npts)     ! r^2 * rho * v
   REAL (KIND=dp) :: dflux_dr(Npts) ! d(r^2 * rho * v)/dr
   REAL (KIND=dp) :: drho_dr(Npts)  ! d(rho)/dr
   REAL (KIND=dp) :: dv_dr(Npts)    ! d(v)/dr
   REAL (KIND=dp) :: grav(Npts)     ! gravitational acceleration dPhi/dr
   REAL (KIND=dp) :: drho_dt(Npts)  ! time derivative of rho
   REAL (KIND=dp) :: dv_dt(Npts)    ! time derivative of v
   INTEGER :: i

   ! Compute the flux r^2 * rho * v for the continuity equation 
   DO i = 1, Npts
      flux(i) = r(i)**2 * rho_in(i) * v_in(i)
   END DO

   ! Gradients 
   dflux_dr = GRADIENT(flux,    r, Npts)
   drho_dr= GRADIENT(rho_in,  r, Npts)
   dv_dr = GRADIENT(v_in,    r, Npts)

   ! Gravitational force 
   IF (use_gravity) THEN
      grav = GRAV_FORCE(rho_in, r, Npts)
   ELSE
      grav = 0.0_dp
   END IF

   ! Time derivatives 
   DO i = 1, Npts
      ! Continuity: d(rho)/dt = -(1/r^2) * d(r^2*rho*v)/dr
      drho_dt(i) = -(1.0_dp / r(i)**2) * dflux_dr(i)

      ! Euler: d(v)/dt = -v*dv/dr - (cs^2/rho)*d(rho)/dr - dPhi/dr
      dv_dt(i) = - v_in(i) * dv_dr(i) - (cs**2 / rho_in(i)) * drho_dr(i) - grav(i)
   END DO

   ! Explicit Euler step 
   rho_out = rho_in + drho_dt * dt
   v_out  = v_in  + dv_dt   * dt

   ! Boundary conditions zero-gradient
   rho_out(1) = rho_out(2)
   v_out(1) = v_out(2)
   rho_out(Npts) = rho_out(Npts-1)
   v_out(Npts) = v_out(Npts-1)

   ! Safety prevent negative densities 
   DO i = 1, Npts
      IF (rho_out(i) <= 0.0_dp) THEN
         WRITE(*,'(A,I5,A,ES12.4)') '  WARNING: Negative density at index ', i, ', rho = ', rho_out(i)
         rho_out(i) = 1.0e-40_dp   ! floor to tiny positive value
      END IF
   END DO

   END SUBROUTINE PRESTEL_STEP_SPHERICAL
!===============================================================================


!===============================================================================
!  Advance the system by one timestep 1D Cartesian

!
!    d(rho)/dt = -d(rho * v)/dx
!    d(v)/dt   = -v * dv/dx - (cs^2/rho) * d(rho)/dx

!===============================================================================
   SUBROUTINE PRESTEL_STEP_CARTESIAN(rho_in, v_in, r, cs, dt, Npts, rho_out, v_out)

   IMPLICIT NONE

   INTEGER,        INTENT(IN)  :: Npts
   REAL (KIND=dp), INTENT(IN)  :: rho_in(Npts), v_in(Npts), r(Npts)
   REAL (KIND=dp), INTENT(IN)  :: cs, dt
   REAL (KIND=dp), INTENT(OUT) :: rho_out(Npts), v_out(Npts)

   REAL (KIND=dp) :: flux(Npts)    ! rho * v
   REAL (KIND=dp) :: dflux_dx(Npts)
   REAL (KIND=dp) :: drho_dx(Npts)
   REAL (KIND=dp) :: dv_dx(Npts)
   REAL (KIND=dp) :: drho_dt(Npts)
   REAL (KIND=dp) :: dv_dt(Npts)
   INTEGER :: i

   ! Compute the flux r^2 * rho * v for the continuity equation 
   ! (In Cartesian mode r is just x and r^2 cancels, but we still call GRADIENT)
   DO i = 1, Npts
      flux(i) = rho_in(i) * v_in(i)
   END DO

   ! Gradients with periodic wrap-around at boundaries 
   ! Interior: standard centred difference
   DO i = 2, Npts-1
      dflux_dx(i) = (flux(i+1)   - flux(i-1))   / (r(i+1) - r(i-1))
      drho_dx(i) = (rho_in(i+1) - rho_in(i-1)) / (r(i+1) - r(i-1))
      dv_dx(i) = (v_in(i+1)   - v_in(i-1))   / (r(i+1) - r(i-1))
   END DO
   ! Boundary: periodic wrap (point 1 neighbours Npts-1 on the left; point Npts neighbours 2 on the right)
   dflux_dx(1)    = (flux(2) - flux(Npts-1))    / (r(2) + (r(Npts) - r(Npts-1)) - r(1))
   drho_dx(1)     = (rho_in(2) - rho_in(Npts-1))  / (r(2) + (r(Npts) - r(Npts-1)) - r(1))
   dv_dx(1)       = (v_in(2) - v_in(Npts-1))    / (r(2) + (r(Npts) - r(Npts-1)) - r(1))
   dflux_dx(Npts) = (flux(2) - flux(Npts-1))    / (r(2) + (r(Npts) - r(Npts-1)) - r(1))
   drho_dx(Npts)  = (rho_in(2) - rho_in(Npts-1))  / (r(2) + (r(Npts) - r(Npts-1)) - r(1))
   dv_dx(Npts)    = (v_in(2) - v_in(Npts-1))    / (r(2) + (r(Npts) - r(Npts-1)) - r(1))

   ! Time derivatives 
   DO i = 1, Npts
      drho_dt(i) = -dflux_dx(i)
      dv_dt(i)   = -v_in(i) * dv_dx(i) - (cs**2 / rho_in(i)) * drho_dx(i)
   END DO

   ! Explicit Euler step 
   rho_out = rho_in + drho_dt * dt
   v_out  = v_in   + dv_dt   * dt

   rho_out(1) = rho_out(Npts-1)
   v_out(1) = v_out(Npts-1)
   rho_out(Npts) = rho_out(2)
   v_out(Npts)   = v_out(2)

   WHERE (rho_out <= 1.0e-35_dp) rho_out = 1.0e-35_dp

   ! Safet floor on density 
   DO i = 1, Npts
      IF (rho_out(i) <= 0.0_dp) rho_out(i) = 1.0e-40_dp
   END DO

   END SUBROUTINE PRESTEL_STEP_CARTESIAN
!===============================================================================


!===========================
   END MODULE PRESTEL
!===========================
