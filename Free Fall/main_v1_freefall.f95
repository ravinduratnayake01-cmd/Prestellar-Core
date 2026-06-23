!==============================================================================

PROGRAM MAIN_PRESTEL

   USE CONSTANTS
   USE PRESTEL

   IMPLICIT NONE

   INTEGER, PARAMETER :: Npts = 200
   REAL(KIND=dp), PARAMETER :: CFL_fac = 0.5_dp

   CALL RUN_TEST_FREEFALL()


CONTAINS

!==============================================================================
   SUBROUTINE RUN_TEST_FREEFALL()

      IMPLICIT NONE

      REAL(KIND=dp), PARAMETER :: cs_ff = 1.0_dp      ! nearly pressureless [cm/s]
      REAL(KIND=dp), PARAMETER :: rho0 = 1.0e-19_dp   ! uniform density [g/cm^3]
      REAL(KIND=dp), PARAMETER :: Rmin = 1.0e13_dp    ! inner radius [cm]
      REAL(KIND=dp), PARAMETER :: Rmax = 1.0e16_dp    ! outer radius [cm]

      REAL(KIND=dp) :: r(Npts), rho(Npts), v(Npts)
      REAL(KIND=dp) :: rho_new(Npts), v_new(Npts)
      REAL(KIND=dp) :: t_ff, dt, dt_cfl, t
      LOGICAL :: cfl_ok
      INTEGER :: i, step, total_steps, imid

      WRITE(*,*) '  Free-fall collapse '

      ! Analytic free-fall time
      t_ff = SQRT(3.0_dp * pi / (32.0_dp * GNewton * rho0))
      WRITE(*,'(A,ES12.4,A)') ' rho0  = ', rho0, ' g/cm^3'
      WRITE(*,'(A,ES12.4,A)') ' t_ff  = ', t_ff, ' s'
      WRITE(*,'(A,ES12.4,A)') ' t_ff  = ', t_ff / yr_to_s, ' yr'

      CALL RADIAL_GRID(Npts, Rmin, Rmax, 'lin', r)

      DO i = 1, Npts
         rho(i) = rho0
         v(i) = 0.0_dp
      END DO

      CALL CFL_CHECK(v, r, cs_ff, Npts, CFL_fac, 1.0e30_dp, dt_cfl, cfl_ok)

      dt = MIN(dt_cfl, t_ff / 1000.0_dp)
      ! Limit dt, can never exceed t_ff/1000, stops CFL giving a timestep larger than the free-fall time
      total_steps = INT(0.3_dp * t_ff / dt)

      WRITE(*,'(A,ES12.4,A)') ' dt_CFL = ', dt, ' s'
      WRITE(*,'(A,ES12.4,A)') '  dt used = ', dt, ' s'
      WRITE(*,'(A,ES12.4,A)') '  t_ff = ', t_ff, ' s'
      WRITE(*,'(A,I8)')  ' Steps  = ', total_steps

      ! Write initial state
      OPEN(10, FILE='freefall_initial.dat', STATUS='REPLACE')
      WRITE(10,'(A)') '# r[cm]   rho[g/cm3]   v[cm/s]'
      DO i = 1, Npts
         WRITE(10,'(3ES18.8)') r(i), rho(i), v(i)
      END DO
      CLOSE(10)

      ! Time loop
      t = 0.0_dp
      DO step = 1, total_steps
         CALL PRESTEL_STEP(rho, v, r, cs_ff, dt, Npts, rho_new, v_new)
         rho = rho_new
         v = v_new
         t = t + dt
      END DO

      ! Write final state
      OPEN(11, FILE='freefall_final.dat', STATUS='REPLACE')
      WRITE(11,'(A)') '# r[cm]   rho[g/cm3]   v[cm/s]'
      DO i = 1, Npts
         WRITE(11,'(3ES18.8)') r(i), rho(i), v(i)
      END DO
      CLOSE(11)

      imid = Npts / 2
      WRITE(*,*)
      WRITE(*,'(A,ES12.4,A)') ' Time reached = ', t / t_ff,  ' * t_ff'
      WRITE(*,'(A,ES12.4,A)') ' Central rho (final) = ', rho(imid), ' g/cm^3'
      WRITE(*,'(A,ES12.4)')   ' Density increase = ', rho(imid) / rho0

   END SUBROUTINE RUN_TEST_FREEFALL

END PROGRAM MAIN_PRESTEL
