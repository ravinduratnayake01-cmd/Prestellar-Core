PROGRAM MAIN_PRESTEL

   USE CONSTANTS
   USE PRESTEL

   IMPLICIT NONE

   ! Common parameters
   INTEGER, PARAMETER :: Npts = 200
   REAL(KIND=dp), PARAMETER :: CFL_fac = 0.4_dp

   CALL RUN_TEST1_ACOUSTIC()
   CALL RUN_TEST2_JEANS()
   CALL RUN_TEST3_FREEFALL()
   CALL RUN_STABILITY_ANALYSIS()

CONTAINS

!  TEST 1 — Acoustic wave propagation
!==============================================================================
   SUBROUTINE RUN_TEST1_ACOUSTIC()

      IMPLICIT NONE

      REAL(KIND=dp), PARAMETER :: T_gas = 10.0_dp        ! Temperature [K]
      REAL(KIND=dp), PARAMETER :: rho0 = 1.0e-20_dp      ! Background density [g/cm^3]
      REAL(KIND=dp), PARAMETER :: A_frac = 0.01_dp       ! Perturbation fraction
      REAL(KIND=dp), PARAMETER :: Rmin = 1.0_dp          ! Grid start [cm]
      REAL(KIND=dp), PARAMETER :: Rmax = 1.0e3_dp        ! Grid end [cm]

      REAL(KIND=dp) :: r(Npts), rho(Npts), v(Npts)
      REAL(KIND=dp) :: rho_new(Npts), v_new(Npts)
      REAL(KIND=dp) :: rho_init(Npts)
      REAL(KIND=dp) :: cs, k, lambda, T_period, dt, t, A
      REAL(KIND=dp) :: dt_cfl, rms_err, norm
      LOGICAL :: cfl_ok
      INTEGER :: i, step, total_steps

      WRITE(*,*) '  TEST 1: Acoustic wave'

      CALL RADIAL_GRID(Npts, Rmin, Rmax, 'lin', r)
      cs = SOUND_SPEED(T_gas)
      A = A_frac * rho0

      lambda = Rmax - Rmin          ! One wavelength fills the domain
      k = twopi / lambda
      T_period = lambda / cs

      ! Initial right-going plane wave
      DO i = 1, Npts
         rho(i) = rho0 + A * COS(k * r(i))
         v(i) = cs * (A / rho0) * COS(k * r(i))
      END DO
      rho_init = rho

      ! CFL timestep
      CALL CFL_CHECK(v, r, cs, Npts, CFL_fac, 1.0e30_dp, dt_cfl, cfl_ok)
      dt = dt_cfl
      total_steps = INT(T_period / dt) + 1

      ! Write initial profile
      OPEN(20, FILE='test1_initial.dat', STATUS='REPLACE')
      WRITE(20,'(A)') '# r[cm]   rho[g/cm3]   v[cm/s]'
      DO i = 1, Npts; WRITE(20,'(3ES18.8)') r(i), rho(i), v(i); END DO
      CLOSE(20)

      ! Time integration — use PRESTEL_STEP but with grav=0 trick
      ! gravity is zero because cs is the only force here,
      ! and GRAV_FORCE_SUB for uniform rho0 << perturbation is negligible
      t = 0.0_dp
      DO step = 1, total_steps
         CALL PRESTEL_STEP(rho, v, r, cs, dt, Npts, rho_new, v_new)
         rho = rho_new
         v   = v_new
         t   = t + dt
      END DO

      ! Write final profile
      OPEN(21, FILE='test1_final.dat', STATUS='REPLACE')
      WRITE(21,'(A)') '# r[cm]   rho[g/cm3]   v[cm/s]'
      DO i = 1, Npts; WRITE(21,'(3ES18.8)') r(i), rho(i), v(i); END DO
      CLOSE(21)

      ! RMS error between final and initial density
      rms_err = 0.0_dp; norm = 0.0_dp
      DO i = 1, Npts
         rms_err = rms_err + (rho(i) - rho_init(i))**2
         norm = norm + rho_init(i)**2
      END DO
      rms_err = SQRT(rms_err / norm)

   END SUBROUTINE RUN_TEST1_ACOUSTIC

!  TEST 2 — Jeans stability: sub-Jeans oscillates, super-Jeans collapses
!==============================================================================
   SUBROUTINE RUN_TEST2_JEANS()

      IMPLICIT NONE

      REAL(KIND=dp), PARAMETER :: T_gas = 10.0_dp
      REAL(KIND=dp), PARAMETER :: rho_sub = 1.0e-21_dp      ! sub-Jeans  [g/cm^3]
      REAL(KIND=dp), PARAMETER :: rho_sup = 1.0e-18_dp      ! super-Jeans [g/cm^3]
      REAL(KIND=dp), PARAMETER :: Rmin = 1.0e14_dp          ! ~0.003 pc [cm]
      REAL(KIND=dp), PARAMETER :: Rmax = 1.0e17_dp          ! ~0.03 pc [cm]
      REAL(KIND=dp), PARAMETER :: A_frac = 0.05_dp
      INTEGER, PARAMETER :: Nsteps = 2000

      REAL(KIND=dp) :: r(Npts), rho(Npts), v(Npts)
      REAL(KIND=dp) :: rho_new(Npts), v_new(Npts)
      REAL(KIND=dp) :: cs, dt, dt_cfl, t
      REAL(KIND=dp) :: lambda_J, rho0_run, M_init, M_final
      LOGICAL :: cfl_ok
      INTEGER :: i, step, icase

      WRITE(*,*) '  TEST 2: Jeans stability (sub, super)'

      cs = SOUND_SPEED(T_gas)

      DO icase = 1, 2
         IF (icase == 1) THEN
            rho0_run = rho_sub
         ELSE
            rho0_run = rho_sup
         END IF

         ! Jeans length
         lambda_J = cs * SQRT(pi / (GNewton * rho0_run))
         WRITE(*,'(A,ES12.4,A)') ' rho0 = ', rho0_run, ' g/cm^3'
         WRITE(*,'(A,ES12.4,A)') ' Jeans length = ', lambda_J/AU_to_cm, ' AU'

         CALL RADIAL_GRID(Npts, Rmin, Rmax, 'log', r)

         ! Uniform density + small radial perturbation
         DO i = 1, Npts
            rho(i) = rho0_run * (1.0_dp + A_frac * COS(pi * r(i) / Rmax))
            v(i)   = 0.0_dp
         END DO

         M_init = TOTAL_MASS(rho, r, Npts)

         CALL CFL_CHECK(v, r, cs, Npts, CFL_fac, 1.0e30_dp, dt_cfl, cfl_ok)
         dt = dt_cfl
         WRITE(*,'(A,ES12.4,A)') ' dt_CFL = ', dt, ' s'

         t = 0.0_dp
         DO step = 1, Nsteps
            CALL PRESTEL_STEP(rho, v, r, cs, dt, Npts, rho_new, v_new)
            rho = rho_new
            v   = v_new
            t   = t + dt
         END DO

         M_final = TOTAL_MASS(rho, r, Npts)
         WRITE(*,'(A,ES12.4,A)') ' Central rho at t_end = ', rho(Npts/2), ' g/cm^3'
         WRITE(*,'(A,ES12.4)')   ' Mass conservation dM/M = ', &
                                  ABS(M_final - M_init) / M_init

         IF (icase == 1) THEN
            OPEN(30, FILE='test2_subjeans.dat', STATUS='REPLACE')
         ELSE
            OPEN(30, FILE='test2_superjeans.dat', STATUS='REPLACE')
         END IF
         WRITE(30,'(A)') '# r[cm]   rho[g/cm3]   v[cm/s]'
         DO i = 1, Npts; WRITE(30,'(3ES18.8)') r(i), rho(i), v(i); END DO
         CLOSE(30)

      END DO


   END SUBROUTINE RUN_TEST2_JEANS


!  TEST 3 — Free fall collapse
!==============================================================================
   SUBROUTINE RUN_TEST3_FREEFALL()

      IMPLICIT NONE

      ! To simulate free-fall, use a very small cs
      REAL(KIND=dp), PARAMETER :: cs_ff = 1.0_dp        ! nearly pressureless [cm/s]
      REAL(KIND=dp), PARAMETER :: rho0 = 1.0e-19_dp
      REAL(KIND=dp), PARAMETER :: Rmin = 1.0e13_dp
      REAL(KIND=dp), PARAMETER :: Rmax = 1.0e16_dp

      REAL(KIND=dp) :: r(Npts), rho(Npts), v(Npts)
      REAL(KIND=dp) :: rho_new(Npts), v_new(Npts)
      REAL(KIND=dp) :: t_ff, dt, dt_cfl, t
      LOGICAL :: cfl_ok
      INTEGER :: i, step, total_steps, imid

      WRITE(*,*) '  TEST 3: Free-fall collapse (near-pressureless)'

      ! Analytic free-fall time: t_ff = sqrt(3*pi / (32*G*rho0))
      t_ff = SQRT(3.0_dp * pi / (32.0_dp * GNewton * rho0))
      WRITE(*,'(A,ES12.4,A)') ' rho0 = ', rho0, ' g/cm^3'
      WRITE(*,'(A,ES12.4,A)') ' t_ff = ', t_ff, ' s'
      WRITE(*,'(A,ES12.4,A)') ' t_ff = ', t_ff / yr_to_s, ' yr'

      CALL RADIAL_GRID(Npts, Rmin, Rmax, 'lin', r)

      ! Uniform sphere at rest
      DO i = 1, Npts
         rho(i) = rho0
         v(i)   = 0.0_dp
      END DO

      CALL CFL_CHECK(v, r, cs_ff, Npts, CFL_fac, 1.0e30_dp, dt_cfl, cfl_ok)
      dt = dt_cfl
      total_steps = INT(0.3_dp * t_ff / dt)   ! integrate to 30% of t_ff

      WRITE(*,'(A,ES12.4,A)') ' dt_CFL = ', dt, ' s'
      WRITE(*,'(A,I8)') ' Steps = ', total_steps

      ! Write initial state
      OPEN(40, FILE='test3_initial.dat', STATUS='REPLACE')
      WRITE(40,'(A)') '# r[cm]   rho[g/cm3]   v[cm/s]'
      DO i = 1, Npts; WRITE(40,'(3ES18.8)') r(i), rho(i), v(i); END DO
      CLOSE(40)

      t = 0.0_dp
      DO step = 1, total_steps
         CALL PRESTEL_STEP(rho, v, r, cs_ff, dt, Npts, rho_new, v_new)
         rho = rho_new
         v = v_new
         t = t + dt
      END DO

      imid = Npts / 2
      WRITE(*,'(A,ES12.4,A)') ' Time reached: ', t / t_ff, ' * t_ff'
      WRITE(*,'(A,ES12.4,A)') ' Central rho at t_end = ', rho(imid), ' g/cm^3'
      WRITE(*,'(A,ES12.4)')   ' Density increase factor = ', rho(imid) / rho0

      IF (rho(imid) > rho0 * 1.5_dp) THEN
         WRITE(*,*) '  PASS: Density increased — collapse proceeding.'
      ELSE
         WRITE(*,*) '  WARN: Density barely changed — check cs_ff or grid.'
      END IF

      OPEN(41, FILE='test3_final.dat', STATUS='REPLACE')
      WRITE(41,'(A)') '# r[cm]   rho[g/cm3]   v[cm/s]'
      DO i = 1, Npts; WRITE(41,'(3ES18.8)') r(i), rho(i), v(i); END DO
      CLOSE(41)

   END SUBROUTINE RUN_TEST3_FREEFALL

END PROGRAM MAIN_PRESTEL
