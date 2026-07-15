      program Arfcc
      implicit none
      integer nmax,n,i,j,k,l,inc,m,bs,rec
      real*8 bohr,sgm,amass
      parameter (nmax=800000)
      parameter(bohr=0.5292d0)
      parameter(sgm=3.4d0/bohr) 
      parameter(amass=40d0*1836d0)
      real*8 hxx,hyy,hzz,tempK
      real*8 x(3*nmax),v(3,nmax)
      real*8 xp(4),yp(4),zp(4),cunit,treg,kb,g(3)
      real*8 kinetic_energy
      cunit=2.0d0**(1d0/6)*sgm*
     &  sqrt(2d0)*0.975d0
      rec=0
      write(*,*)cunit

! boxsize
      bs=15
! boxsize        
      hxx=bs*cunit*bohr
      hyy=bs*cunit*bohr
      hzz=bs*cunit*bohr

 
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

      
      inc=0
      do i=0,bs-1
      do j=0,bs-1
      do k=0,bs-1
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
      write(*,*)n

! 初速度      
      treg=100.d0
      kb=1.d0/11605.d0/27.2116d0
      do i=1,n
         call gauss(g)
         v(:,i)=sqrt(kb*treg/amass)*g(:)
      enddo
      kinetic_energy=0.d0
      do i=1,n
       kinetic_energy=kinetic_energy+0.5d0*amass*sum(v(:,i)**2)
      enddo
       tempK=kinetic_energy*2.d0/(3.d0*n)*11605.d0*27.2116d0
      write(*,*)tempK     
!

      open(10,file='init.dat')
      write(10,*)n
      do i=1,n
        write(10,'(a,1x,i5,6e15.7)') 'Ar',i,
     &       x(3*i-2)*bohr,x(3*i-1)*bohr,x(3*i)*bohr,
     &       v(1,i),v(2,i),v(3,i)
      enddo
      write(10,'(3e24.15)')hxx,0d0,0d0
      write(10,'(3e24.15)')0d0,hyy,0d0
      write(10,'(3e24.15)')0d0,0d0,hzz
      close(10)

      open(11,file='first.xyz')

          write(11,*)n
          write(11,'(a,3(3e15.7),a,a)')
     &       'Lattice="',hxx,0.0,0.0,0.0,hyy,0.0,0.0,0.0,hzz,'" ',
     &       'Properties=species:S:1:id:I:1:pos:R:3:tempK:R:1'
           do m=1,n
            tempK=amass/2*(v(1,m)**2+v(2,m)**2+v(3,m)**2)
     &        *27.2116*11605d0
            write(11,'(a2,1x,i5,4e15.7)') 'Ar',m,
     &       x(3*m-2)*bohr,x(3*m-1)*bohr,x(3*m)*bohr,tempK
          enddo
          close(11)
      end


 ! maxwell分布
      subroutine gauss(g)
      implicit none
      real*8 g(3),u1,u2
      integer i

      do i=1,3
       call random_number(u1)
       call random_number(u2)
       if(u1<1d-12) u1=1d-12
       g(i)=sqrt(-2.d0*log(u1))*cos(2.d0*3.14159265d0*u2)
      enddo
      return
      end
 !



