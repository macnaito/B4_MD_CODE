      program Arfcc
      implicit none
      integer nmax,n,i,j,k,l,inc,m
      real*8 bohr,sgm
      parameter (nmax=30000)
      parameter(bohr=0.5292d0)
      parameter(sgm=3.4d0/bohr) 

      real*8 hxx,hyy,hzz
      real*8 x(3*nmax),v(3*nmax)
      real*8 xp(4),yp(4),zp(4),cunit

      hxx=100d0
      hyy=100d0
      hzz=100d0
 
      xp(1)=0d0
      yp(1)=0d0
      zp(1)=0d0
      xp(2)=0.5d0
      yp(2)=0.5d0
      zp(2)=0d0
      xp(3)=0.5d0
      yp(3)=0d0
      zp(3)=0.5d0
      xp(4)=0d0
      yp(4)=0.5d0
      zp(4)=0.5d0

      cunit=2.0d0**(1d0/6)*sgm*sqrt(2d0)*0.98
      write(*,*) cunit
      inc=0
      do i=0,5
      do j=0,5
      do k=0,5
        do l=1,4
          inc=inc+1
          x(3*inc-2)=(xp(l)+dble(i))*cunit
          x(3*inc-1)=(yp(l)+dble(j))*cunit
          x(3*inc  )=(zp(l)+dble(k))*cunit
        enddo
      enddo
      enddo
      enddo
      n=inc
      do i=1,n
        v(3*i-2)= 0d0
        v(3*i-1)= 0d0
        v(3*i  )= 0d0
      enddo
      

      
      open(10,file='init.dat')
      write(10,*)n
      do i=1,n
        write(10,'(a,i5,6e15.7)') 'Ar',i,
     &       x(3*i-2)*bohr,x(3*i-1)*bohr,x(3*i)*bohr,
     &       v(3*i-2),v(3*i-1),v(3*i)
      enddo
      write(10,'(3e24.15)')hxx,0d0,0d0
      write(10,'(3e24.15)')0d0,hyy,0d0
      write(10,'(3e24.15)')0d0,0d0,hzz

      close(10)

!      open(11,file='init.xyz')
!          write(11,*)n
!          write(11,'(a,3(3e15.7),a,a)')
!     &       'Lattice="',hxx,0.0,0.0,0.0,hyy,0.0,0.0,0.0,hzz,'" ',
!     &       'Properties=species:S:1:id:I:1:pos:R:3:tempK:R:1'
!         do m=1,n
!            write(11,'(a2,i5,4e15.7)')'Ar',m,
!     &       x(3*m-2)*bohr,x(3*m-1)*bohr,x(3*m)*bohr,0.0d0
!          enddo
!      close(11)

      end






