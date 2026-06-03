!   << 3-dimension >>
!
!   Dynamics of n atoms interacting with the Lennard-Jones potential.
      program md_lj_pbc
      implicit none
      integer nmax,rec
      parameter(nmax=5000)
      integer i,j,k,m,n,maxstep,itemp
      real*8 x(3*nmax),v(3*nmax),dfdx(3*nmax)
      real*8 f,ekin,dt
      real*8 hxx,hyy,hzz,temp,dummy,tempK
      real*8 amass,bohr,hx,hy,hz,Tk
      character*2 lsp(nmax)
      character*40 filename2
      real*8 Q,tau,Treg,xi,gkbt
      parameter(amass=40d0*1836d0)
      parameter(bohr=0.5292d0)
! maxstep
      parameter(maxstep=2000)
      real*8 T(maxstep)
      real*8 poten(maxstep)
      real*8 xii(maxstep)

!ファイルの読み込み
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
       read(10,*)hx,dummy,temp
       read(10,*)temp,hy,temp
       read(10,*)temp,temp,hz
       hxx=hx/bohr
       hyy=hy/bohr
       hzz=hz/bohr
      close(10)
      ekin=0d0
      do j=1,3*n
          ekin=ekin+0.5d0*amass*v(j)**2
      enddo

      filename2='out000.xyz'
      k=0
      Tk=0.d0
      rec=0
      xi=0.d0
! time_step
      dt=41d0
! 目標温度
      Treg=75d0
      tau=95d0*dt
      gkbt=(3.d0*dble(n))*Treg/(27.2116*11605d0)
      Q=gkbt*tau**2

      call pot(f,dfdx,x,n,hxx,hyy,hzz)

      do i=1,maxstep
        xi=xi+0.5d0*dt*(2.d0*ekin-gkbt)/Q
        ekin=0d0
        do j=1,3*n
          v(j)=(v(j)-(dt*0.5d0)*dfdx(j)/amass)
     &         /(1+dt*0.5d0*xi)
          ekin=ekin+0.5d0*amass*v(j)**2 
        enddo

! 周期境界条件
        do j=1,n
          x(3*j-2)=x(3*j-2)+dt*v(3*j-2)
          if (x(3*j-2).gt.hxx) then
            x(3*j-2)=x(3*j-2)-hxx
          else if (x(3*j-2).lt.0d0) then
            x(3*j-2)=x(3*j-2)+hxx
          endif 
          x(3*j-1)=x(3*j-1)+dt*v(3*j-1)
          if (x(3*j-1).gt.hyy) then
            x(3*j-1)=x(3*j-1)-hyy
          else if (x(3*j-1).lt.0d0) then
            x(3*j-1)=x(3*j-1)+hyy
          endif
          x(3*j  )=x(3*j  )+dt*v(3*j  )
          if (x(3*j  ).gt.hzz) then
            x(3*j  )=x(3*j  )-hzz
          else if (x(3*j  ).lt.0d0) then
            x(3*j  )=x(3*j  )+hzz
          endif
        enddo  

        call pot(f,dfdx,x,n,hxx,hyy,hzz)
        poten(i)=f

        xi=xi+0.5d0*dt*(2.d0*ekin-gkbt)/Q

        ekin=0d0
        do j=1,3*n
         v(j)=(v(j)-dt*0.5d0*dfdx(j)/amass)
     &     /(1+dt*0.5d0*xi)
         ekin=ekin+0.5d0*amass*v(j)**2 
        enddo


        Tk=ekin*2d0/(3d0*dble(n))*
     &    27.2116*11605d0
        T(i)=Tk
        xii(i)=xi
        write(*,*) Tk
        


!-------Note: 1 atomic unit of time = 2.42d-17 sec
!        write(*,*)dt*i*2.42d-17,ekin,f,totalenergy

        if(mod(i,100).eq.0)then
          k=k+1
          write(filename2(4:6),'(i3.3)')k

          open(11,file=filename2)
          write(11,*)n
          write(11,'(a,3(3e15.7),a,a)')
     &       'Lattice="',hx,0.0,0.0,0.0,hy,0.0,0.0,0.0,hz,'" ',
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

      if (rec==1) then
       open (10,file='final50x.dat')
        do i=1,3*n
         write(10,*) x(i)
        enddo
       close(10)
       open (10,file='final50v.dat')
        do i=1,3*n
         write(10,*) v(i)
        enddo
       close(10)
      endif

      open (12,file='t.dat')
      do i=1,maxstep
        write(12,*) T(i),xii(i)*1000000,poten(i)
      enddo
      close(12)
    

          

      

      endprogram md_lj_pbc

      subroutine pot(f,dfdx,x,n,hxx,hyy,hzz)
      implicit none

      integer n,i,j
      real*8 f,x(3*n),dfdx(3*n)
      real*8 xij,yij,zij,r2,factor
      real*8 sgm,eps,sgm12,sgm6
      real*8 hxx,hyy,hzz,cutoff
      parameter(sgm=3.4d0/0.5292d0)
      parameter(eps=120d0/11605d0/27.2116d0)
      parameter(sgm12=sgm**12,sgm6=sgm**6)
      parameter(cutoff=2.5d0*sgm)
      f=0d0
c-----potential
      do i=1,n-1
        do j=i+1,n
          xij=x(3*i-2)-x(3*j-2)
          yij=x(3*i-1)-x(3*j-1)
          zij=x(3*i  )-x(3*j  )
          xij=xij-hxx*dnint(xij/hxx)
          yij=yij-hyy*dnint(yij/hyy)
          zij=zij-hzz*dnint(zij/hzz)
          r2=xij**2+yij**2+zij**2
          if (r2 < cutoff**2) then
            f=f+4d0*eps*(sgm12/r2**6-sgm6/r2**3)
          endif
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
            xij=xij-hxx*dnint(xij/hxx)
            yij=yij-hyy*dnint(yij/hyy)
            zij=zij-hzz*dnint(zij/hzz)
            r2=xij**2+yij**2+zij**2
            if (r2 < cutoff**2) then
              factor=4d0*eps*
     &        (-12d0*sgm12/r2**7+6d0*sgm6/r2**4)
              dfdx(3*i-2)=dfdx(3*i-2)+factor*xij
              dfdx(3*i-1)=dfdx(3*i-1)+factor*yij
              dfdx(3*i  )=dfdx(3*i  )+factor*zij
            endif
          endif
        enddo
      enddo

      return
      end