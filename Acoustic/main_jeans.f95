! Test 2: Jeans instability
! Dispersion relation: omega^2 = cs^2 * k^2 - 4*pi*G*rho0

! a) omega^2 < 0  (lambda > lambdaJ): exponential perturbation
! b) omega^2 > 0  (lambda < lambdaJ): oscillation
! c) omega^2 = 0  (lambda = lambdaJ): critical wavelength

!===============================================================================

PROGRAM JEANS_TEST

   USE CONSTANTS
   USE PRESTEL
   USE JEANS_ANALYSIS

   IMPLICIT NONE

   INTEGER,       PARAMETER :: Npts = 400
   INTEGER,       PARAMETER :: n_waves = 2
   REAL(KIND=dp), PARAMETER :: Rmin = 0.0_dp
   REAL(KIND=dp), PARAMETER :: T_gas = 10.0_dp
   REAL(KIND=dp), PARAMETER :: rho0 = 1.0e-20_dp
   REAL(KIND=dp), PARAMETER :: A_frac = 1.0e-3_dp
   REAL(KIND=dp), PARAMETER :: CFL_fac = 0.4_dp
   INTEGER,       PARAMETER :: N_snap = 6             ! density snapshots

   REAL(KIND=dp), DIMENSION(Npts) :: x, rho, v, rho_new, v_new, drho1_dx
   REAL(KIND=dp) :: cs, kJ, lJ, A
   REAL(KIND=dp) :: lam, k_t, Rmax, dr
   REAL(KIND=dp) :: sigma, omega2, omega, gfac
   REAL(KIND=dp) :: dt, t, rho1_max
   INTEGER :: i, step, nsteps, nsnap

   cs = SOUND_SPEED(T_gas)
   kJ = JEANS_WAVENUMBER(rho0, cs)
   lJ = JEANS_LENGTH(rho0, cs)
   A  = A_frac * rho0

   WRITE(*,*)
   WRITE(*,*) 'Test 2: Jeans instability'
   WRITE(*,'(A,ES12.4,A)') '  cs = ', cs, ' cm/s'
   WRITE(*,'(A,ES12.4,A)') '  rho0 = ', rho0, ' g/cm^3'
   WRITE(*,'(A,ES12.4,A)') '  lambdaJ = ', lJ/parsec_to_cm, ' pc'
   WRITE(*,'(A,ES12.4,A)') '  tJ  = ', JEANS_TIME(rho0)/yr_to_s, ' yr'
   WRITE(*,*)



! a): omega^2 < 0 — exponential growth
!===============================================================================
   WRITE(*,*) 'a): omega^2 < 0, lambda = 2*lambdaJ '

   lam = 2.0_dp * lJ
   k_t = twopi / lam
   sigma = GROWTH_RATE(k_t, rho0, cs)
   gfac = 4.0_dp * pi * GNewton / k_t**2

   Rmax = n_waves * lam
   x = RADIAL_GRID(Npts, Rmin, Rmax, 'lin')
   dr = x(2) - x(1)
   dt = CFL_fac * dr / cs

   ! growing mode: rho1 = A*cos(kx), v1 = (sigma/rho0/k)*A*sin(kx)
   DO i = 1, Npts
      rho(i) = rho0 + A * COS(k_t * x(i))
      v(i) = (sigma / (rho0 * k_t)) * A * SIN(k_t * x(i))
   END DO

   ! run for 2 e-folding times
   nsteps = INT(2.0_dp / (sigma * dt)) + 1
   nsnap  = MAX(1, nsteps / N_snap)

   WRITE(*,'(A,ES12.4,A)') '  sigma  = ', sigma, ' s^-1'
   WRITE(*,'(A,ES12.4,A)') '  1/sigma  = ', 1.0_dp/sigma/yr_to_s, ' yr'
   WRITE(*,'(A,ES12.4,A)') '  dt  = ', dt/yr_to_s, ' yr'
   WRITE(*,'(A,I9)')        '  nsteps  = ', nsteps
   WRITE(*,*)

   OPEN(20, FILE='jeans_growth_a.dat', STATUS='REPLACE')
   WRITE(20,'(A)') '# t[yr]  max(rho1)_numerical  max(rho1)_analytic=A*exp(sigma*t)'

   t = 0.0_dp
   DO step = 1, nsteps
      rho1_max = MAXVAL(ABS(rho - rho0))
      WRITE(20,'(3ES16.6)') t/yr_to_s, rho1_max, A*EXP(sigma*t)

      IF (MOD(step, nsnap) == 0) THEN
         WRITE(21,'(A,ES12.4,A)') '# t = ', t/yr_to_s, ' yr'
         DO i = 1, Npts; WRITE(21,'(2ES16.6)') x(i), rho(i); END DO
         WRITE(21,*)
      END IF

      CALL PRESTEL_STEP_CARTESIAN(rho, v, x, cs, dt, Npts, rho_new, v_new)
      DO i = 2, Npts-1
         drho1_dx(i) = (rho(i+1) - rho(i-1)) / (x(i+1) - x(i-1))
      END DO
      drho1_dx(1)    = (rho(2) - rho(Npts-1)) / (x(2)-x(1) + x(Npts)-x(Npts-1))
      drho1_dx(Npts) = drho1_dx(1)
      v_new = v_new + gfac * drho1_dx * dt

      rho = rho_new;  v = v_new;  t = t + dt
   END DO

   CLOSE(20);  CLOSE(21)

   rho1_max = MAXVAL(ABS(rho - rho0))
   WRITE(*,*) 'Case a):'
   WRITE(*,'(A,ES12.4)') '  final max(rho1) numerical = ', rho1_max
   WRITE(*,'(A,ES12.4)') '  final max(rho1) analytic  = ', A*EXP(sigma*t)
   WRITE(*,'(A,F8.4)')   '  ratio   = ', rho1_max / (A*EXP(sigma*t))
   WRITE(*,*)


! Case b): omega^2 > 0 — oscillation
!===============================================================================
   WRITE(*,*) ' Case (b): omega^2 > 0, lambda = 0.5*lambdaJ '

   lam = 0.5_dp * lJ
   k_t = twopi / lam
   omega2 = DISPERSION_OMEGA2(k_t, rho0, cs)
   omega = SQRT(omega2)
   gfac= 4.0_dp * pi * GNewton / k_t**2

   Rmax = n_waves * lam
   x = RADIAL_GRID(Npts, Rmin, Rmax, 'lin')
   dr = x(2) - x(1)
   dt = CFL_fac * dr / cs

   ! right-going wave
   DO i = 1, Npts
      rho(i) = rho0 + A * COS(k_t * x(i))
      v(i) = (omega / (rho0 * k_t)) * A * COS(k_t * x(i))
   END DO

   ! run for one full oscillation period T = 2pi/omega
   nsteps = INT(twopi / (omega * dt)) + 1
   nsnap = MAX(1, nsteps / N_snap)

   WRITE(*,'(A,ES12.4,A)') '  omega  = ', omega, ' rad/s'
   WRITE(*,'(A,ES12.4,A)') '  T = 2pi/omega= ', twopi/omega/yr_to_s, ' yr'
   WRITE(*,'(A,ES12.4,A)') '  dt = ', dt/yr_to_s, ' yr'
   WRITE(*,'(A,I9)')       '  nsteps = ', nsteps
   WRITE(*,*)

   OPEN(22, FILE='jeans_growth_b.dat', STATUS='REPLACE')
   WRITE(22,'(A)') '# t[yr]  max(rho1)_numerical  initial_A  (amplitude should stay ~A)'

   t = 0.0_dp
   DO step = 1, nsteps
      rho1_max = MAXVAL(ABS(rho - rho0))
      WRITE(22,'(3ES16.6)') t/yr_to_s, rho1_max, A

      IF (MOD(step, nsnap) == 0) THEN
         WRITE(23,'(A,ES12.4,A)') '# t = ', t/yr_to_s, ' yr'
         DO i = 1, Npts; WRITE(23,'(2ES16.6)') x(i), rho(i); END DO
         WRITE(23,*)
      END IF

      CALL PRESTEL_STEP_CARTESIAN(rho, v, x, cs, dt, Npts, rho_new, v_new)
      DO i = 2, Npts-1
         drho1_dx(i) = (rho(i+1) - rho(i-1)) / (x(i+1) - x(i-1))
      END DO
      drho1_dx(1) = (rho(2) - rho(Npts-1)) / (x(2)-x(1) + x(Npts)-x(Npts-1))
      drho1_dx(Npts) = drho1_dx(1)
      v_new = v_new + gfac * drho1_dx * dt

      rho = rho_new;  v = v_new;  t = t + dt
   END DO

   CLOSE(22);  CLOSE(23)

   rho1_max = MAXVAL(ABS(rho - rho0))
   WRITE(*,*) 'Case b):'
   WRITE(*,'(A,ES12.4)') '  amplitude at t=0 = ', A
   WRITE(*,'(A,ES12.4)') '  amplitude at t=T = ', rho1_max
   WRITE(*,'(A,F8.4)')   '  ratio (want ~1)  = ', rho1_max / A
   WRITE(*,*)

! Case c): omega^2 = 0 — critical wavelength lambdaJ
!===============================================================================
   WRITE(*,*) ' Case (c): omega^2 = 0, lambda = lambdaJ (critical) '

   lam = lJ
   k_t = twopi / lam
   gfac = 4.0_dp * pi * GNewton / k_t**2

   Rmax = n_waves * lam
   x = RADIAL_GRID(Npts, Rmin, Rmax, 'lin')
   dr = x(2) - x(1)
   dt = CFL_fac * dr / cs

   DO i = 1, Npts
      rho(i) = rho0 + A * COS(k_t * x(i))
      v(i)   = 0.0_dp
   END DO

   ! run for one Jeans time
   nsteps = INT(JEANS_TIME(rho0) / dt)

   WRITE(*,'(A,ES12.4,A)') '  lambdaJ      = ', lJ/parsec_to_cm, ' pc'
   WRITE(*,'(A,ES12.4)')   '  omega^2      = ', DISPERSION_OMEGA2(k_t, rho0, cs)
   WRITE(*,'(A,ES12.4,A)') '  tJ           = ', JEANS_TIME(rho0)/yr_to_s, ' yr'
   WRITE(*,'(A,ES12.4,A)') '  dt           = ', dt/yr_to_s, ' yr'
   WRITE(*,'(A,I9)')       '  nsteps       = ', nsteps
   WRITE(*,*)

   OPEN(24, FILE='jeans_marginal.dat', STATUS='REPLACE')
   WRITE(24,'(A)') '# t[yr]  max(rho1)  initial_A  (should remain constant)'

   t = 0.0_dp
   DO step = 1, nsteps
      rho1_max = MAXVAL(ABS(rho - rho0))
      WRITE(24,'(3ES16.6)') t/yr_to_s, rho1_max, A

      CALL PRESTEL_STEP_CARTESIAN(rho, v, x, cs, dt, Npts, rho_new, v_new)
      DO i = 2, Npts-1
         drho1_dx(i) = (rho(i+1) - rho(i-1)) / (x(i+1) - x(i-1))
      END DO
      drho1_dx(1) = (rho(2) - rho(Npts-1)) / (x(2)-x(1) + x(Npts)-x(Npts-1))
      drho1_dx(Npts) = drho1_dx(1)
      v_new = v_new + gfac * drho1_dx * dt

      rho = rho_new;  v = v_new;  t = t + dt
   END DO

   CLOSE(24)

   rho1_max = MAXVAL(ABS(rho - rho0))
   WRITE(*,*) 'Case c):'
   WRITE(*,'(A,ES12.4)') '  initial amplitude = ', A
   WRITE(*,'(A,ES12.4)') '  final   amplitude = ', rho1_max
   WRITE(*,'(A,F10.4)')  '  drift / A         = ', (rho1_max - A) / A

END PROGRAM JEANS_TEST
