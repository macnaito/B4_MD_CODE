      program Arfcc
      implicit none
      integer nmax,n,i,j,k,l,inc,m,bs,rec,aa
      real*8 bohr,sgm,amass,ekin,Tk,a
      parameter (nmax=30000)
      parameter(bohr=0.5292d0)
      parameter(sgm=3.4d0/bohr) 
      parameter(amass=40d0*1836d0)
      real*8 hxx,hyy,hzz,tempK,sgmv,sgmvv
      real*8 x(3*nmax),v(3*nmax)
      real*8 xp(4),yp(4),zp(4),cunit

      cunit=2.0d0**(1d0/6)*sgm*
     &  sqrt(2d0)*0.975d0
      rec=0

! boxsize
      bs=10
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

! 初速度      
      do i=1,n
        v(3*i-2)=0d0
        v(3*i-1)=0d0
        v(3*i  )=0d0
      enddo
      rec=0
      if (rec==1) then
       open(10,file='final50x.dat')
       do i=1,3*n
            read(10,*)x(i)
       enddo
       close(10)
       open(11,file='final50v.dat')
       do i=1,3*n
            read(11,*)v(i)
      enddo
      close(11)

!格子欠陥
      do j=1,51     
       if (j==1) exit
       call random_number(a)
       aa=int(a*n)+1
       write(*,*)aa
       do i=aa,n-1
            x(3*i-2)=x(3*(i+1)-2)
            x(3*i-1)=x(3*(i+1)-1)
            x(3*i  )=x(3*(i+1)  )
            v(3*i-2)=v(3*(i+1)-2)
            v(3*i-1)=v(3*(i+1)-1)
            v(3*i  )=v(3*(i+1)  )
       enddo
       n=n-1
      enddo

!重心速度をひく     
      sgmv=0d0
      sgmvv=0d0
      do i=1,3*n
            sgmv=sgmv+v(i)
      enddo
      write(*,*)sgmv
      sgmv=sgmv/(3d0*n)
      do i=1,3*n
            v(i)=v(i)-sgmv
      enddo
      do i=1,3*n
            sgmvv=sgmvv+v(i)
      enddo
      do i=1,3*n
        v(i)=(50d0/50.52d0)**0.5*v(i)
      enddo
      write(*,*) sgmvv
      end if
!重心速度

      ekin=0d0
      do j=1,3*n
          ekin=ekin+0.5d0*amass*v(j)**2
      enddo
      Tk=ekin*2d0/(3d0*dble(n))*
     &    27.2116*11605d0
      write(*,*) Tk


      
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

      open(11,file='first.xyz')

          write(11,*)n
          write(11,'(a,3(3e15.7),a,a)')
     &       'Lattice="',hxx,0.0,0.0,0.0,hyy,0.0,0.0,0.0,hzz,'" ',
     &       'Properties=species:S:1:id:I:1:pos:R:3:tempK:R:1'
           do m=1,n
            tempK=amass/2*(v(3*m-2)**2+v(3*m-1)**2+v(3*m)**2)
     &        *27.2116*11605d0
            write(11,'(a2,i5,4e15.7)') 'Ar',m,
     &       x(3*m-2)*bohr,x(3*m-1)*bohr,x(3*m)*bohr,tempK
          enddo
          close(11)
      end



