      program pbc_lvBAOAB_parrinello
! PBC lc lengevine Parrinello-Rahman(isotropic)
      implicit none
      integer nmax
      parameter(nmax=50000)
      integer i,j,k,m,n,maxstep,itemp
      real*8 x(3*nmax),v(3*nmax),dfdx(3*nmax)
      real*8 f,ekin,dt
      real*8 hxx,hyy,hzz,temp,dummy,tempK
      real*8 amass,bohr,hx,hy,hz,TK
      character*2 lsp(nmax)
      character*40 filename2
      parameter(amass=40d0*1836d0)
      parameter(bohr=0.5292d0)
      real*8 L,vL,aL,W,Lold,scale
      real*8 virial,p,kb,pext
! time_step      
      parameter(maxstep=2000)
      real*8 T(maxstep)
      real*8 gpa(maxstep)
! rang
      real*8 gamma,Treg
      real*8 sigma
      real*8 gauss
      external gauss
      dt=41d0*5
      gamma=1.d-4
      kb=1.d0/(27.2116*11605d0)

!ファイルの読みこみ      
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
!ファイルの読み込み

      filename2='out000.xyz'
      k=0
      pext=1.01325d-4
      L=hx
      vL=0.d0
      W=1.d8


!計算start    
      call pot (f,dfdx,x,n,hx,hy,hz,virial)
      ekin=0d0
        do j=1,n
          ekin=ekin+0.5d0*amass*
     &       (v(3*j-2)**2+v(3*j-1)**2+v(3*j)**2)
        enddo
      Tk=ekin*2d0/(3d0*dble(n))*27.2116*11605d0
      p=(n*kb*Tk+virial/3.d0)/(hx*hy*hz)*29421.d0
      aL=(p-pext)*L/W
      vL=vL+0.50*aL*dt  
!計算start 
      do i=1,maxstep
        Treg=50d0     
!        if (i<(Treg-50d00)*200)then
!         Treg=50d0+0.005*i
!        endif

        sigma=sqrt(kb*Treg/amass*(1d0-exp(-gamma*dt)**2))

        do j=1,3*n
          v(j)=v(j)+(dt/2d0)*(dfdx(j)/amass)
        enddo

!圧力の計算
        
        Lold=L
        L=L+vL*0.5d0*dt
        scale=L/Lold
        hx=hx*scale
        hy=hy*scale
        hz=hz*scale
        do j=1,3*n
          x(j)=x(j)*scale
          v(j)=v(j)/scale
        enddo
!位置更新
        do j=1,n
          x(3*j-2)=x(3*j-2)+0.5d0*dt*v(3*j-2)
          x(3*j-1)=x(3*j-1)+0.5d0*dt*v(3*j-1)
          x(3*j  )=x(3*j  )+0.5d0*dt*v(3*j  )
        enddo  
!

        do j=1,3*n
          v(j)=exp(-gamma*dt)*v(j)+sigma*gauss()
        enddo

        Lold=L
        L=L+vL*0.5d0*dt
        scale=L/Lold
        hx=hx*scale
        hy=hy*scale
        hz=hz*scale
        do j=1,3*n
          x(j)=x(j)*scale
          v(j)=v(j)/scale
        enddo

!位置更新
        do j=1,n
          x(3*j-2)=x(3*j-2)+0.5d0*dt*v(3*j-2)          
          x(3*j-1)=x(3*j-1)+0.5d0*dt*v(3*j-1)
          x(3*j  )=x(3*j  )+0.5d0*dt*v(3*j  )
        enddo  
!

!周期境界条件
        do j=1,n
          if (x(3*j-2).gt.hx) then
            x(3*j-2)=x(3*j-2)-hx
          else if (x(3*j-2).lt.0d0) then
            x(3*j-2)=x(3*j-2)+hx
          endif 
          if (x(3*j-1).gt.hy) then
            x(3*j-1)=x(3*j-1)-hy
          else if (x(3*j-1).lt.0d0) then
            x(3*j-1)=x(3*j-1)+hy
          endif
          if (x(3*j  ).gt.hz) then
            x(3*j  )=x(3*j  )-hz
          else if (x(3*j  ).lt.0d0) then
            x(3*j  )=x(3*j  )+hz
          endif
        enddo
!周期境界条件

        call pot(f,dfdx,x,n,hx,hy,hz,virial)
        do j=1,3*n
          v(j)=v(j)+(dt/2d0)*(dfdx(j)/amass)
        enddo


!温度圧力計算        
        ekin=0d0
        do j=1,n
          ekin=ekin+0.5d0*amass*
     &       (v(3*j-2)**2+v(3*j-1)**2+v(3*j)**2)
        enddo
        Tk=ekin*2d0/(3d0*dble(n))*
     &    27.2116*11605d0
        T(i)=Tk
        p=(n*kb*Tk+virial/3.d0)/(hx*hy*hz)*29421.d0
        gpa(i)=p
        aL=(p-pext)*L/W
        vL=vL+0.5d0*aL*dt
!温度圧力計算        

!ファイルの書き出し
        if(mod(i,100).eq.0)then
          k=k+1
          write(filename2(4:6),'(i3.3)')k

          open(11,file=filename2)
          write(11,*)n
          write(11,'(a,3(3e15.7),a,a)')
     &       'Lattice="',hx*bohr,0.0,0.0,0.0,hy*bohr,
     &           0.0,0.0,0.0,hz*bohr,'" ',
     &       'Properties=species:S:1:id:I:1:pos:R:3:tempK:R:1'
           do m=1,n
            tempK=amass/2*(v(3*m-2)**2+v(3*m-1)**2+v(3*m)**2)
     &        *27.2116*11605d0
            write(11,'(a2,i5,4e15.7)')lsp(m),m,
     &       x(3*m-2)*bohr,x(3*m-1)*bohr,x(3*m)*bohr,tempK
          enddo
          close(11)
          TK=ekin*315775.0d0/(1.5d0*n)
        endif
!ここまで
      write(*,*) Tk,p

      enddo    
 

!温度のグラフ      
      open (12,file='t.dat')
      do i=1,maxstep
        write(12,*) T(i),gpa(i)
      enddo
      close(12)
!ここまで
      open(10,file='lv50_pr1d-4_8788.dat')
      write(10,*)n
      do i=1,n
       write(10,'(a,i5,6e15.7)') 'Ar',i,
     &       x(3*i-2)*bohr,x(3*i-1)*bohr,x(3*i)*bohr,
     &       v(3*i-2),v(3*i-1),v(3*i)
      enddo
      write(10,'(3e24.15)')hx*bohr,0d0,0d0
      write(10,'(3e24.15)')0d0,hx*bohr,0d0
      write(10,'(3e24.15)')0d0,0d0,hz*bohr
      close(10)

      endprogram 

!linked cell
      subroutine pot(f,dfdx,x,n,hx,hy,hz,virial)
      implicit real*8(a-h,o-z)
      integer mx,my,mz,hx_lc,hy_lc,hz_lc
      real*8 x(3*n),dfdx(3*n),f,factor
      integer hyz_lc,hxyz_lc,i
      integer ishiftx,ishifty,ishiftz
      parameter(sgm=3.4d0/0.5292d0)
      parameter(eps=120d0/11605d0/27.2116d0)
      parameter(sgm12=sgm**12,sgm6=sgm**6)
      parameter(cutoff=2.5d0*sgm)
      integer,allocatable,dimension(:):: lshd,lscl
      real*8 virial

      virial=0.d0
      do i=1,3*n
        dfdx(i)=0d0
      enddo

      hx_lc=max(int(hx/cutoff),1) !x座標のセル数
      hy_lc=max(int(hy/cutoff),1)
      hz_lc=max(int(hz/cutoff),1)
      hyz_lc=hy_lc*hz_lc     !セルのyz平面の数
      hxyz_lc=hyz_lc*hx_lc     !セルの総数
      hx_cell=hx/hx_lc     !各方向のセルの長さ
      hy_cell=hy/hy_lc
      hz_cell=hz/hz_lc
      allocate(lshd(hxyz_lc),lscl(n))
      lshd=0
      do i=1,n
        if(x(3*i-2).lt.0d0 .or. x(3*i-2).ge.hx .or.
     &     x(3*i-1).lt.0d0 .or. x(3*i-1).ge.hy .or.
     &     x(3*i  ).lt.0d0 .or. x(3*i  ).ge.hz )then
          write(*,*)'Error! i, x,y,z=',
     &     i,x(3*i-2),x(3*i-1),x(3*i)
          stop
        endif
        mx=int(x(3*i-2)/hx_cell)
        my=int(x(3*i-1)/hy_cell)
        mz=int(x(3*i  )/hz_cell)
        mx=min(max(mx,0),hx_lc-1)
        my=min(max(my,0),hy_lc-1)
        mz=min(max(mz,0),hz_lc-1)
        m=mx*hyz_lc+my*hz_lc+mz+1
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
        m=mx*hyz_lc+my*hz_lc+mz+1 !3つのdoloopで全てのセルを走査
        if (lshd(m)==0)cycle !空のセルはスキップ
        do kuz=-kuzmax,kuzmax
        do kuy=-kuymax,kuymax
        do kux=-kuxmax,kuxmax
          m1x=mx+kux
          call get_ishift(hx_lc,m1x,ishiftx)
          m1y=my+kuy
          call get_ishift(hy_lc,m1y,ishifty)
          m1z=mz+kuz
          call get_ishift(hz_lc,m1z,ishiftz)
          m1=m1x*hyz_lc+m1y*hz_lc+m1z+1 !走査するセル番号
          if (lshd(m1)==0) cycle !走査するセルがからならスキップ
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
                  dfdx(3*i-2)=dfdx(3*i-2)-factor*rijx
                  dfdx(3*i-1)=dfdx(3*i-1)-factor*rijy
                  dfdx(3*i  )=dfdx(3*i  )-factor*rijz
                  dfdx(3*j-2)=dfdx(3*j-2)+factor*rijx
                  dfdx(3*j-1)=dfdx(3*j-1)+factor*rijy
                  dfdx(3*j  )=dfdx(3*j  )+factor*rijz
                  virial=virial-factor*rij2
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

!ガウス乱数 box muller 平均０分散１
      function gauss()
      implicit none
      real*8 gauss
      real*8 u1,u2

      call random_number(u1)
      call random_number(u2)
      gauss=sqrt(-2.d0*log(u1))
     & *cos(2.d0*3.14159265d0*u2) 
      
      return
      end