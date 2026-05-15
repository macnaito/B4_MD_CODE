!   << 1-dimension >>
!
!   Dynamics of 3 particles with mass, amass=1, interacting with
!    the harmonic potential.
!
      program md_3particles
      implicit none

      integer, parameter :: nmax=10
      integer :: n, maxstep
      integer :: i, j, k
      real(8) :: dt, amass
      real(8) :: ekin, totalenergy, f

      real(8) :: x0(nmax), v0(nmax)
      real(8) :: x(nmax), v(nmax), dfdx(nmax)

      n = 3
      maxstep = 5000
      dt = 0.01d0
      amass = 1.0d0

      x0(1) = -2.0d0
      v0(1) =  0.5d0
      x0(2) =  0.0d0
      v0(2) =  0.1d0
      x0(3) =  2.0d0
      v0(3) = -0.6d0

      do i = 1, n
        x(i) = x0(i)
        v(i) = v0(i)
      end do

      call pot(f, dfdx, x, n)

      do i = 1, maxstep

        do j = 1, n
          v(j) = v(j) + (dt/2.0d0) * ( -dfdx(j) / amass )
        end do

        do j = 1, n
          x(j) = x(j) + dt * v(j)
        end do

        call pot(f, dfdx, x, n)

        do j = 1, n
          v(j) = v(j) + (dt/2.0d0) * ( -dfdx(j) / amass )
        end do

        ekin = 0.0d0
        do j = 1, n
          ekin = ekin + 0.5d0 * amass * v(j)**2
        end do

        totalenergy = f + ekin

        write(*,*) dt*i, (x(k), k=1,n), totalenergy

      end do

      end program md_3particles
      
!---------------------------------------------------------------
      subroutine pot(f, dfdx, x, n)
      implicit none

      integer, intent(in) :: n
      real(8), intent(in) :: x(n)
      real(8), intent(out) :: dfdx(n)
      real(8), intent(out) :: f

      integer :: i

      f = 0.0d0

      do i = 1, n-1
        f = f + (x(i+1) - x(i) - 1.0d0)**2
      end do

      do i = 1, n
        if (i == 1) then
          dfdx(i) = 2.0d0 * (x(i) - x(i+1) + 1.0d0)
        elseif (i == n) then
          dfdx(i) = 2.0d0 * (x(i) - x(i-1) - 1.0d0)
        else
          dfdx(i) = 2.0d0 * (x(i) - x(i-1) - 1.0d0)
     &              + 2.0d0 * (x(i) - x(i+1) + 1.0d0)
        end if
      end do

      return
      end