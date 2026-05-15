      implicit real*8(a-h,o-z)
      parameter(nmax=3000)
      real*8 xp(4),yp(4),zp(4)
      real*8 x(3*nmax),v(3*nmax)
!-----sigmal of LJ pot. in atomic unit
!-----Note: 1 atomic unit of length = 0.5292Angstrom
      parameter(bohr=0.5292d0)
      parameter(sgm=3.4d0/bohr)

!-----MD box
      hxx=200d0
      hyy=200d0
      hzz=200d0

!-----face-centered cubic: primitive unit that contains 4 atoms
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

      cunit=2d0**(1d0/6)*sgm*sqrt(2d0)*0.98

!-----base plane
      inc=0
      do i=-5,5
      do j=-5,5
      do k=-5,-1
        do l=1,4
          inc=inc+1
          x(3*inc-2)=(xp(l)+dble(i))*cunit
          x(3*inc-1)=(yp(l)+dble(j))*cunit
          x(3*inc  )=(zp(l)+dble(k))*cunit
        enddo
      enddo
      enddo
      enddo
      n1=inc

!-----small cube
      do i=-4,-1
      do j=-1,2
      do k= 0,3
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

      if(n.ge.nmax)then
        write(*,*)'Too large n=',n
        stop
      endif

!-----shift all atoms
      do i=1,n
        x(3*i-2)=x(3*i-2)+hxx/2
        x(3*i-1)=x(3*i-1)+hyy/2
        x(3*i  )=x(3*i )+hzz/2
      enddo

!-----set velocities to the base plane
      do i=1,n1
        v(3*i-2)= 0d0
        v(3*i-1)= 0d0
        v(3*i  )= 0d0
      enddo

!-----set vx=100m/sec to the small cube
!-----Note: 1 atomic unit of velocity = 2.19d6 m/sec
      do i=n1+1,n
        v(3*i-2)= 100d0/2.19d6
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

      end
