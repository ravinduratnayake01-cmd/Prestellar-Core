! Test 1: Acoustic wave in 1D Cartesian geometry
! No gravity. Uniform background rho0, small perturbation A*cos(kx).

!===============================================================================
PROGRAM ACOUSTIC_TEST

   USE CONSTANTS
   USE PRESTEL

   IMPLICIT NONE

   !  parameters 
   INTEGER,       PARAMETER :: Npts  = 200
   INTEGER,       PARAMETER :: n_waves = 4
   REAL(KIND=dp), PARAMETER :: T_gas = 10.0_dp           ! K
   REAL(KIND=dp), PARAMETER :: rho0 = 1.0e-20_dp         ! g/cm^3
   REAL(KIND=dp), PARAMETER :: A_frac = 0.01_dp          ! perturbation amplitude / rho0
   REAL(KIND=dp), PARAMETER :: lambda = 1.0e15_dp        ! cm  (~0.3 AU)
   REAL(KIND=dp), PARAMETER :: Rmin = 0.0_dp
   REAL(KIND=dp), PARAMETER :: Rmax  = n_waves * lambda  ! Make the domain size match a reasonable wavelength
   REAL(KIND=dp), PARAMETER :: CFL_fac = 0.4_dp

   !  local 
   REAL(KIND=dp), DIMENSION(Npts) :: x, rho, v, rho_new, v_new, rho_init
   REAL(KIND=dp) :: cs, k, T_period, dt, t, A
   REAL(KIND=dp) :: rms, norm
   INTEGER :: i, step, nsteps, nhalf

   ! setup
   x = RADIAL_GRID(Npts, Rmin, Rmax, 'lin')
   cs = SOUND_SPEED(T_gas)
   A = A_frac * rho0
   k = twopi / lambda
   T_period = lambda / cs
   dt = CFL_fac * (x(2) - x(1)) / cs
   nsteps = INT(T_period / dt)
   nhalf  = nsteps / 2

   WRITE(*,*)
   WRITE(*,*) 'Test 1: Acoustic wave'
   WRITE(*,'(A,ES12.4,A)') '  cs  = ', cs,  ' cm/s'
   WRITE(*,'(A,ES12.4,A)') '  lambda = ', lambda, ' cm'
   WRITE(*,'(A,ES12.4,A)') '  T_period = ', T_period/yr_to_s,' yr'
   WRITE(*,'(A,ES12.4,A)') '  dt = ', dt/yr_to_s, ' yr'
   WRITE(*,'(A,I8)')       '  nsteps = ', nsteps
   WRITE(*,'(A,F8.1)')     '  pts/wave  = ', lambda / (x(2)-x(1))
   WRITE(*,*)

   ! initial conditions: right-going wave
   ! rho1 = A*cos(kx),  v1 = cs*(A/rho0)*cos(kx)
   DO i = 1, Npts
      rho(i) = rho0 + A * COS(k * x(i))
      v(i)   = cs * (A / rho0) * COS(k * x(i))
   END DO
   rho_init = rho

   ! write t=0
   OPEN(10, FILE='wave_t0.dat', STATUS='REPLACE')
   WRITE(10,'(A)') '# x[cm]  rho[g/cm3]  v[cm/s]'
   DO i = 1, Npts
      WRITE(10,'(3ES16.6)') x(i), rho(i), v(i)
   END DO
   CLOSE(10)

   ! time loop
   t = 0.0_dp
   DO step = 1, nsteps
      IF (step == nhalf) THEN
         OPEN(11, FILE='wave_thalf.dat', STATUS='REPLACE')
         WRITE(11,'(A)') '# x[cm]  rho[g/cm3]  v[cm/s]'
         DO i = 1, Npts
            WRITE(11,'(3ES16.6)') x(i), rho(i), v(i)
         END DO
         CLOSE(11)
      END IF
      CALL PRESTEL_STEP_CARTESIAN(rho, v, x, cs, dt, Npts, rho_new, v_new)
      rho = rho_new
      v = v_new
      t = t + dt
   END DO

   ! write t~T
   OPEN(12, FILE='wave_tT.dat', STATUS='REPLACE')
   WRITE(12,'(A)') '# x[cm]  rho[g/cm3]  v[cm/s]'
   DO i = 1, Npts
      WRITE(12,'(3ES16.6)') x(i), rho(i), v(i)
   END DO
   CLOSE(12)

   ! analytic solution at t
   OPEN(13, FILE='wave_analytic.dat', STATUS='REPLACE')
   WRITE(13,'(A)') '# x[cm]  rho_analytic[g/cm3]  v_analytic[cm/s]  at t=T'
   DO i = 1, Npts
      WRITE(13,'(3ES16.6)') x(i), &
         rho0 + A * COS(k*x(i) - cs*k*t), &
         cs * (A/rho0) * COS(k*x(i) - cs*k*t)
   END DO
   CLOSE(13)

   ! error at t~T vs t=0
   rms  = 0.0_dp
   norm = 0.0_dp
   DO i = 1, Npts
      rms  = rms  + (rho(i) - rho_init(i))**2
      norm = norm + rho_init(i)**2
   END DO
   rms = SQRT(rms / norm)

   WRITE(*,*) 'Results:'
   WRITE(*,'(A,ES12.4)') '  RMS(rho_T - rho_0) / RMS(rho_0) = ', rms
   WRITE(*,*) '  (small => wave shape preserved over one period)'
   WRITE(*,*)
   WRITE(*,*) 'Output files: wave_t0.dat  wave_thalf.dat  wave_tT.dat  wave_analytic.dat'

END PROGRAM ACOUSTIC_TEST
