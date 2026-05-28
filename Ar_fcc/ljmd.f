!   << 3-dimension >>
!
!   Dynamics of n atoms interacting with the Lennard-Jones potential.
!
      program md_lj
      implicit none

      integer nmax
      parameter(nmax=2000)

      integer i,j,k,m,n,maxstep,itemp
      real*8 x(3*nmax),v(3*nmax),dfdx(3*nmax)
      real*8 f,ekin,totalenergy,dt
      real*8 hxx,hyy,hzz,temp,dummy,tempK
      real*8 amass,bohr

      character*2 lsp(nmax)
      character*40 filename,filename2

!-----mass of the atom in atomic unit
      parameter(amass=40d0*1836d0)
!-----Note: 1 atomic unit of length = 0.5292Angstrom
      parameter(bohr=0.5292d0)

!-----time step in atomic unit.
      dt=41d0*5d0
!-----total number of steps
      maxstep=2000

!-----read the initial atomic positions and velocities
      open(10,file='init.dat')
      read(10,*)n
      do i=1,n
        read(10,*)lsp(i),itemp,
     &  x(3*i-2),x(3*i-1),x(3*i),
     &  v(3*i-2),v(3*i-1),v(3*i)

        x(3*i-2)=x(3*i-2)/bohr
        x(3*i-1)=x(3*i-1)/bohr
        x(3*i  )=x(3*i  )/bohr
      enddo
      read(10,*)hxx,dummy,temp
      read(10,*)temp,hyy,temp
      read(10,*)temp,temp,hzz
      close(10)

      filename2='out000.xyz'

      k=0

!-----velocity-Verlet algorithm
      call pot(f,dfdx,x,n)

      do i=1,maxstep

        do j=1,3*n
          v(j)=v(j)+(dt/2d0)*(-dfdx(j)/amass)
        enddo

        do j=1,3*n
          x(j)=x(j)+dt*v(j)
        enddo  

        call pot(f,dfdx,x,n)

        do j=1,3*n
          v(j)=v(j)+(dt/2d0)*(-dfdx(j)/amass)
        enddo

        ekin=0d0
        do j=1,3*n
          ekin=ekin+0.5d0*amass*v(j)**2
        enddo

        totalenergy=f+ekin

!-------Note: 1 atomic unit of time = 2.42d-17 sec
        write(*,*)dt*i*2.42d-17,ekin,f,totalenergy

        if(mod(i,100).eq.0)then

          k=k+1
          write(filename(4:6),'(i3.3)')k
          write(filename2(4:6),'(i3.3)')k

c-----chem3d
!          open(10,file=filename)
!          write(10,*)n
!          do m=1,n
!            write(10,'(a2,i5,6e15.7)')lsp(m),m,
!     &      x(3*m-2)*bohr,x(3*m-1)*bohr,x(3*m)*bohr,
!     &      v(3*m-2),v(3*m-1),v(3*m)
!          enddo
!          write(10,'(3e24.15)')hxx,0d0,0d0
!          write(10,'(3e24.15)')0d0,hyy,0d0
!          write(10,'(3e24.15)')0d0,0d0,hzz
!          close(10)

c-----xyz
          open(11,file=filename2)
          write(11,*)n
          write(11,'(a,3(3e15.7),a,a)')
     &       'Lattice="',hxx,0.0,0.0,0.0,hyy,0.0,0.0,0.0,hzz,'" ',
     &       'Properties=species:S:1:id:I:1:pos:R:3:tempK:R:1'
           do m=1,n
            tempK=amass/2*(v(3*m-2)**2+v(3*m-1)**2+v(3*m)**2)
     &        *27.2116*11605d0
            write(11,'(a2,i5,4e15.7)')lsp(m),m,
     &       x(3*m-2)*bohr,x(3*m-1)*bohr,x(3*m)*bohr,tempK
          enddo
          close(11)

        endif

      enddo

      end

      subroutine pot(f,dfdx,x,n)
      implicit none

      integer n,i,j
      real*8 f,x(3*n),dfdx(3*n)
      real*8 xij,yij,zij,r2,factor
      real*8 sgm,eps,sgm12,sgm6
!-----atomic unit is used
      parameter(sgm=3.4d0/0.5292d0)
      parameter(eps=120d0/11605d0/27.2116d0)
      parameter(sgm12=sgm**12,sgm6=sgm**6)

      f=0d0

c-----potential
      do i=1,n-1
        do j=i+1,n
          xij=x(3*i-2)-x(3*j-2)
          yij=x(3*i-1)-x(3*j-1)
          zij=x(3*i  )-x(3*j  )
          r2=xij**2+yij**2+zij**2
          f=f+4d0*eps*(sgm12/r2**6-sgm6/r2**3)
        enddo
      enddo

c-----force
      do i=1,n
        dfdx(3*i-2)=0d0
        dfdx(3*i-1)=0d0
        dfdx(3*i  )=0d0

        do j=1,n
          if(j.ne.i)then
            xij=x(3*i-2)-x(3*j-2)
            yij=x(3*i-1)-x(3*j-1)
            zij=x(3*i  )-x(3*j  )
            r2=xij**2+yij**2+zij**2

            factor=4d0*eps*
     &      (-12d0*sgm12/r2**7+6d0*sgm6/r2**4)

            dfdx(3*i-2)=dfdx(3*i-2)+factor*xij
            dfdx(3*i-1)=dfdx(3*i-1)+factor*yij
            dfdx(3*i  )=dfdx(3*i  )+factor*zij
          endif
        enddo
      enddo

      return
      end