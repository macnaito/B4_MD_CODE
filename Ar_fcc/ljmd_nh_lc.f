!   << 3-dimension >>
!
!   Dynamics of n atoms interacting with the Lennard-Jones potential.
      program md_lj_pbc
      implicit none
      integer nmax,rec
      parameter(nmax=30000)
      integer i,j,k,m,n,maxstep,itemp
      real*8 x(3*nmax),v(3*nmax),dfdx(3*nmax)
      real*8 f,ekin,dt
      real*8 hxx,hyy,hzz,temp,dummy,tempK
      real*8 amass,bohr,hx,hy,hz,Tk
      character*2 lsp(nmax)
      character*40 filename2
      real*8 Q,tau,xi,gkbt,hamil,eta
      parameter(amass=40d0*1836d0)
      parameter(bohr=0.5292d0)
! maxstep
      parameter(maxstep=10000)
      real*8 T(maxstep)
      real*8 H(maxstep)
      real*8 Treg(maxstep)

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
       read(10,*)hxx,dummy,temp
       read(10,*)temp,hyy,temp
       read(10,*)temp,temp,hzz
       hx=hxx/bohr
       hy=hyy/bohr
       hz=hzz/bohr
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
      eta=0.d0
! time_step
      dt=41d0*5

      call pot(f,dfdx,x,n,hx,hy,hz)
      do i=1,maxstep
! 目標温度
!       Treg=100d0
        Treg(i)=50d0+300d0*dble(i)/maxstep
        tau=40d0*dt
        gkbt=(3.d0*dble(n))*Treg(i)/(27.2116*11605d0)
        Q=gkbt*tau**2
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
          if (x(3*j-2).gt.hx) then
            x(3*j-2)=x(3*j-2)-hx
          else if (x(3*j-2).lt.0d0) then
            x(3*j-2)=x(3*j-2)+hx
          endif 
          x(3*j-1)=x(3*j-1)+dt*v(3*j-1)
          if (x(3*j-1).gt.hy) then
            x(3*j-1)=x(3*j-1)-hy
          else if (x(3*j-1).lt.0d0) then
            x(3*j-1)=x(3*j-1)+hy
          endif
          x(3*j  )=x(3*j  )+dt*v(3*j  )
          if (x(3*j  ).gt.hz) then
            x(3*j  )=x(3*j  )-hz
          else if (x(3*j  ).lt.0d0) then
            x(3*j  )=x(3*j  )+hz
          endif
        enddo  
!周期境界条件

! call pot
        call pot(f,dfdx,x,n,hx,hy,hz)

        xi=xi+0.5d0*dt*(2.d0*ekin-gkbt)/Q

        ekin=0d0
        do j=1,3*n
         v(j)=(v(j)-dt*0.5d0*dfdx(j)/amass)
     &     /(1+dt*0.5d0*xi)
         ekin=ekin+0.5d0*amass*v(j)**2 
        enddo

        eta=eta+xi*dt
        hamil=ekin+f+0.5d0*Q*xi**2+gkbt*eta
        Tk=ekin*2d0/(3d0*dble(n))*
     &    27.2116*11605d0
        T(i)=Tk
        H(i)=hamil
        write(*,*) Tk

!-------Note: 1 atomic unit of time = 2.42d-17 sec
!        write(*,*)dt*i*2.42d-17,ekin,f,totalenergy
!ファイルへの書き出し
        if(mod(i,100).eq.0)then
          k=k+1
          write(filename2(4:6),'(i3.3)')k        
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
!ファイルへの書き出し

!記録  
      rec=0
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
        write(12,*) T(i),H(i)
      enddo
      close(12)
!記録

      endprogram md_lj_pbc

      subroutine pot(f,dfdx,x,n,hx,hy,hz)
      implicit real*8(a-h,o-z)
      integer mx,my,mz,hx_lc,hy_lc,hz_lc
      real*8 x(3*n),dfdx(3*n),f,factor
      integer lcyz,lcxyz,i
      integer ishiftx,ishifty,ishiftz
      parameter(sgm=3.4d0/0.5292d0)
      parameter(eps=120d0/11605d0/27.2116d0)
      parameter(sgm12=sgm**12,sgm6=sgm**6)
      parameter(cutoff=2.5d0*sgm)
      integer,allocatable,dimension(:):: lshd,lscl

      do i=1,3*n
        dfdx(i)=0d0
      enddo

      hx_lc=max(int(hx/cutoff),1) !x座標のセル数
      hy_lc=max(int(hy/cutoff),1)
      hz_lc=max(int(hz/cutoff),1)
      lcyz=hy_lc*hz_lc     !セルのyz平面の数
      lcxyz=lcyz*hx_lc     !セルの総数
      hx_cell=hx/hx_lc     !各方向のセルの大きさ
      hy_cell=hy/hy_lc
      hz_cell=hz/hz_lc
      allocate(lshd(lcxyz),lscl(n))
      lshd=0
      do i=1,n
!        if(x(3*i-2).lt.0d0 .or. x(3*i-2).ge.hx .or.
!     &     x(3*i-1).lt.0d0 .or. x(3*i-1).ge.hy .or.
!     &     x(3*i  ).lt.0d0 .or. x(3*i  ).ge.hz )then
!          write(*,*)'Error! i, x,y,z=',
!     &     i,x(3*i-2),x(3*i-1),x(3*i)
!          stop
!        endif
        mx=int(x(3*i-2)/hx_cell)
        my=int(x(3*i-1)/hy_cell)
        mz=int(x(3*i  )/hz_cell)
        mx=min(max(mx,0),hx_lc-1)
        my=min(max(my,0),hy_lc-1)
        mz=min(max(mz,0),hz_lc-1)
        m=mx*lcyz+my*hz_lc+mz+1
        lscl(i)=lshd(m)
        lshd(m)=i
      enddo

      kuxmax=int(cutoff/hx+1)
      kuymax=int(cutoff/hy+1)
      kuzmax=int(cutoff/hz+1)
      f=0d0
      
      do mz=0,hz_lc-1
      do my=0,hy_lc-1
      do mx=0,hx_lc-1
        m=mx*lcyz+my*hz_lc+mz+1
        if (lshd(m)==0)cycle
        do kuz=-kuzmax,kuzmax
        do kuy=-kuymax,kuymax
        do kux=-kuxmax,kuxmax
          m1x=mx+kux
          call get_ishift(hx_lc,m1x,ishiftx)
          m1y=my+kuy
          call get_ishift(hy_lc,m1y,ishifty)
          m1z=mz+kuz
          call get_ishift(hz_lc,m1z,ishiftz)
          m1=m1x*lcyz+m1y*hz_lc+m1z+1
          if (lshd(m1)==0) cycle
          i=lshd(m)
          do while(i>0)
            j=lshd(m1)
            do while(j>0)
              if (i<j) then
                rijx=x(3*i-2)-(x(3*j-2)+ishiftx*hx)
                rijy=x(3*i-1)-(x(3*j-1)+ishifty*hy)
                rijz=x(3*i  )-(x(3*j  )+ishiftz*hz)
                rij2=rijx**2+rijy**2+rijz**2
                if (rij2<cutoff**2) then
                  f=f+4d0*eps*(sgm12/rij2**6-sgm6/rij2**3)
                  factor=4d0*eps*
     &            (-12d0*sgm12/rij2**7+6d0*sgm6/rij2**4)
                  dfdx(3*i-2)=dfdx(3*i-2)+factor*rijx
                  dfdx(3*i-1)=dfdx(3*i-1)+factor*rijy
                  dfdx(3*i  )=dfdx(3*i  )+factor*rijz
                  dfdx(3*j-2)=dfdx(3*j-2)-factor*rijx
                  dfdx(3*j-1)=dfdx(3*j-1)-factor*rijy
                  dfdx(3*j  )=dfdx(3*j  )-factor*rijz
                end if
              endif
              j=lscl(j)
            enddo
            i=lscl(i)
          enddo
        enddo
        enddo
        enddo
      enddo
      enddo
      enddo
      deallocate(lshd,lscl)
      end      


!subrutine ishift
      subroutine get_ishift(lc,i,ishift)
      implicit none
      integer lc,i,ishift
      ishift=0
      if(i.gt.lc-1)then
        do while(i > lc-1)
          ishift=ishift+1
          i=i-lc
        enddo
      elseif(i.lt.0)then
        do while(i < 0)
          ishift=ishift-1
          i=i+lc
        enddo
      else
      endif
      return
      end
!end