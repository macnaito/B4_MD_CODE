      program test_lj
      implicit none

      integer i
      real*8 r,f,dfdr
      real*8 f1,dfdr1

      write(*,*)' r   f  dfdr(true) dfdr(num)'

      do i=40,500
        r=i*0.02d0
        call func(f,dfdr,r)
        call func(f1,dfdr1,r+1d-5)
        write(*,*)r,f,dfdr,(f1-f)/1d-5
      enddo

      end
c----------------------------------------------
      subroutine func(f,dfdr,r)
      implicit none

      real*8 f,dfdr,r
      real*8 eps,sgm

      eps=1d0
      sgm=1d0

      f=4d0*eps*((sgm/r)**12-(sgm/r)**6)

      dfdr=4d0*eps*( -12d0*sgm**12*r**(-13)
     &               + 6d0*sgm**6*r**(-7) )

      return
      end