      program init_config
      implicit none

      integer nmax
      parameter (nmax=10)

      integer n,i
      real*8 bohr,sgm
      real*8 hxx,hyy,hzz
      real*8 x(3*nmax),v(3*nmax)

      parameter (bohr=0.5292d0)
      parameter (sgm=3.4d0/bohr)

      hxx=100d0
      hyy=100d0
      hzz=100d0

      n=3

c----- Atom 1
      x(3*1-2)=-0.7d0*sgm+hxx/2d0
      x(3*1-1)= 0d0*sgm+hyy/2d0
      x(3*1  )= 0d0*sgm+hzz/2d0
      v(3*1-2)=0d0
      v(3*1-1)=0d0
      v(3*1  )=0d0

c----- Atom 2
      x(3*2-2)= 0.7d0*sgm+hxx/2d0
      x(3*2-1)= 0d0*sgm+hyy/2d0
      x(3*2  )= 0d0*sgm+hzz/2d0
      v(3*2-2)=0d0
      v(3*2-1)=0d0
      v(3*2  )=0d0

c----- Atom 3
      x(3*3-2)= 0d0*sgm+hxx/2d0
      x(3*3-1)= 1.2d0*sgm+hyy/2d0
      x(3*3  )= 0d0*sgm+hzz/2d0
      v(3*3-2)=0d0
      v(3*3-1)=0d0
      v(3*3  )=0d0

      open(10,file='init.dat')

      write(10,*) n

      do i=1,n
      write(10,'(a,i5,6e15.7)') 'Ar',i,
     & x(3*i-2)*bohr,x(3*i-1)*bohr,x(3*i)*bohr,
     & v(3*i-2),v(3*i-1),v(3*i)
      enddo
 
      write(10,'(3e24.15)') hxx,0d0,0d0
      write(10,'(3e24.15)') 0d0,hyy,0d0
      write(10,'(3e24.15)') 0d0,0d0,hzz

      close(10)

      end