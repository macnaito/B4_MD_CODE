      program md_lj
! linked_cell pbc parrinelo-rahman 
      implicit none
      integer nmax
      parameter(nmax=5000)
      integer i,j,k,m,n,maxstep,itemp
      real*8 x(3,nmax),v(3,nmax),dfdx(3,nmax)
      real*8 f,dt,kb,Tk,ekin
      real*8 temp,dummys
      real*8 amass,bohr,hx,hy,hz
      character*2 lsp(nmax)
      character*40 filename2
      parameter(amass=40d0*1836d0)
      parameter(bohr=0.5292d0)
! time_step      
      parameter(maxstep=100)
      real*8 T(maxstep)
      

! PR parameter
      real*8 h(3,3),vol,h_inver(3,3)
      real*8 p(3,3),p_reg(3,3),s(3,nmax)
      real*8 vs(3,nmax),as(3,nmax)
      real*8 gh(3,3),ah(3,3),W,sx,sy,sz

!ファイルの読みこみ      
      open(10,file='init.dat')
      read(10,*)n
      do i=1,n
        read(10,*)lsp(i),itemp,
     &  x(1,i),x(2,i),x(3,i),
     &  v(1,i),v(2,i),v(3,i)

        x(1,i)=x(1,i)/bohr
        x(2,i)=x(2,i)/bohr
        x(3,i)=x(3,i)/bohr
      enddo
      read(10,*)hx,dummys,temp
      read(10,*)temp,hy,temp
      read(10,*)temp,temp,hz
      hx=hx/bohr
      hy=hy/bohr
      hz=hz/bohr
      close(10)
!ファイルの読み込み

      h=0d0
      h(1,1)=hx
      h(2,2)=hy
      h(3,3)=hz
      W=1.d8
      kb=1.d0/(27.2116*11605d0)

      call inverse_mass(h,h_inver,vol)
      call pot(f,dfdx,x,n,h,p)
      call pres(v,n,p,vol)

      filename2='out000.xyz'
      k=0
      p_reg(1,1)=0.001
      p_reg(2,2)=0.001
      p_reg(3,3)=0.001
      dt=41d0*5d0

!計算start      
      call pot(f,dfdx,x,n,h,p)
     



      call inverse_mass(h,h_inver,vol)
      ah=matmul((p-p_reg),transpose(h_inver))
      ah=vol*ah/W
      
      do i=1,n
        as(:,i)=matmul(h_inver,dfdx(:,i))/amass
        as(:,i)=as(:,i)-matmul(h_inver,matmul(gh,vs(:,i)))
      enddo


      do i=1,maxstep

        gh=gh+0.5d0*dt*ah
        do j=1,n
          vs(:,j)=vs(:,j)+0.5d0*dt*as(:,j)
        enddo
        h=h+dt*gh
        do j=1,n
          s(:,j)=s(:,j)+dt*vs(:,j)
        enddo

        do j=1,n
         do m=1,3
         if (s(m,j).gt.1.d0) then
          s(m,j)=s(m,j)-1.d0
         elseif (s(m,j).lt.0.d0) then
           s(m,j)=s(m,j)+1.d0
         endif 
         enddo
        enddo

        do j=1,n
          x(:,j)=matmul(h,s(:,j))
        enddo



        call inverse_mass(h,h_inver,vol)
        call pot(f,dfdx,x,n,h,p)
        ekin=0d0
        do j=1,n
         v(:,j)=matmul(h,vs(:,j))
     &        +matmul(gh,s(:,j))
        enddo
        do j=1,n
         ekin=ekin+0.5d0*amass*v(1,j)**2
         ekin=ekin+0.5d0*amass*v(2,j)**2
         ekin=ekin+0.5d0*amass*v(3,j)**2
        enddo
        Tk=ekin*2d0/(3d0*dble(n))*
     &    27.2116*11605d0
        do j=1,3
        do m=1,3
         if(j==m) then
          p(j,m)=(n*kb*Tk+p(j,m))/vol
         else
          p(j,m)=p(j,m)/vol
         endif
        enddo
        enddo


  


!ファイルの書き出し
        if(mod(i,1).eq.0)then
          k=k+1
          write(filename2(4:6),'(i3.3)')k

          open(11,file=filename2)
          write(11,*)n
          write(11,'(a,3(3e15.7),a,a)')
     &       'Lattice="',h(1,1)*bohr,h(2,1)*bohr,h(1,3)*bohr,
     &                   h(2,1)*bohr,h(2,2)*bohr,h(2,3)*bohr,
     &                   h(3,1)*bohr,h(3,2)*bohr,h(3,3)*bohr,'" ',
     &       'Properties=species:S:1:id:I:1:pos:R:3:tempK:R:1'
           do m=1,n
!            tempK=amass/2*(v(3*m-2)**2+v(3*m-1)**2+v(3*m)**2)
!     &        *27.2116*11605d0
            write(11,'(a2,i5,4e15.7)')lsp(m),m,
     &       x(1,m)*bohr,x(2,m)*bohr,x(3,m)*bohr,0d0
          enddo
          close(11)
        endif
!ここまで       

      enddo


!温度のグラフ      
      open (12,file='t.dat')
      do i=1,maxstep
        write(12,*) T(i)
      enddo
      close(12)
!ここまで

      endprogram

! call pot
      subroutine pot(f,dfdx,x,n,h,p)
      implicit none
      integer mx,my,mz,hx_lc,hy_lc,hz_lc,n
      real*8 x(3,n),dfdx(3,n),f,factor
      integer hyz_lc,hxyz_lc,i
      integer ishiftx,ishifty,ishiftz,m
      integer m1x,m1y,m1z,m1,j
      real*8 rijx,rijy,rijz,rij2
      integer kuxmax,kuymax,kuzmax,kux,kuy,kuz
      real*8 sgm,sgm6,sgm12,eps,cutoff
      parameter(sgm=3.4d0/0.5292d0)
      parameter(eps=120d0/11605d0/27.2116d0)
      parameter(sgm12=sgm**12,sgm6=sgm**6)
      parameter(cutoff=2.5d0*sgm)
      integer,allocatable,dimension(:):: lshd,lscl

      real*8 h(3,3)
      real*8 p(3,3)
      real*8 hx_ln,hy_ln,hz_ln    
      real*8 h_inver(3,3)
      real*8 vol,sx,sy,sz

      call inverse_mass(h,h_inver,vol) 

      dfdx=0d0
      p=0d0
      f=0d0
      hx_ln=sqrt(h(1,1)**2+h(2,1)**2+h(3,1)**2)
      hy_ln=sqrt(h(1,2)**2+h(2,2)**2+h(3,2)**2)
      hz_ln=sqrt(h(1,3)**2+h(2,3)**2+h(3,3)**2) 
      hx_lc=max(int(hx_ln/cutoff),1) !x座標のセル数
      hy_lc=max(int(hy_ln/cutoff),1)
      hz_lc=max(int(hz_ln/cutoff),1)
      hyz_lc=hy_lc*hz_lc     !セルのyz平面の数
      hxyz_lc=hyz_lc*hx_lc     !セルの総数
      hx_cell=hx_ln/hx_lc     !各方向のセルの長さ
      hy_cell=hy_ln/hy_lc
      hz_cell=hz_ln/hz_lc
      allocate(lshd(hxyz_lc),lscl(n))
      lshd=0

      do i=1,n
  !正規化
        sx=h_inver(1,1)*x(1,i)+h_inver(1,2)*x(2,i)+
     &   h_inver(1,3)*x(3,i)
        sy=h_inver(2,1)*x(1,i)+h_inver(2,2)*x(2,i)+
     &   h_inver(2,3)*x(3,i)
        sz=h_inver(3,1)*x(1,i)+h_inver(3,2)*x(2,i)+
     &   h_inver(3,3)*x(3,i)  
        mx=int(sx*hx_lc)
        my=int(sy*hy_lc)
        mz=int(sz*hz_lc)
        mx=min(max(mx,0),hx_lc-1)
        my=min(max(my,0),hy_lc-1)
        mz=min(max(mz,0),hz_lc-1)
        m=mx*hyz_lc+my*hz_lc+mz+1
        if(m < 1 .or. m > hxyz_lc) then
         write(*,*) 'BAD M1'
         write(*,*) 'm1=',m
         write(*,*) 'm1x,m1y,m1z=',m1x,m1y,m1z
         write(*,*) 'mx,my,mz=',mx,my,mz
         write(*,*) 'kux,kuy,kuz=',kux,kuy,kuz
         write(*,*) 'hx_lc,hy_lc,hz_lc=',hx_lc,hy_lc,hz_lc
         stop
        endif
        lscl(i)=lshd(m)
        lshd(m)=i
      enddo

      kuxmax=int(cutoff/hx_ln+1)
      kuymax=int(cutoff/hy_ln+1)
      kuzmax=int(cutoff/hz_ln+1)
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
                rijx=x(1,i)-(x(1,j)+
     &           h(1,1)*ishiftx+h(1,2)*ishifty+h(1,3)*ishiftz)
                rijy=x(2,i)-(x(2,j)+
     &           h(2,1)*ishiftx+h(2,2)*ishifty+h(2,3)*ishiftz)
                rijz=x(3,i)-(x(3,j)+
     &           h(3,1)*ishiftx+h(3,2)*ishifty+h(3,3)*ishiftz)
                rij2=rijx**2+rijy**2+rijz**2
                if (rij2<cutoff**2) then
    !ポテンシャル              
                  f=f+4d0*eps*(sgm12/rij2**6-sgm6/rij2**3)
                  factor=24d0*eps*
     &            (-2d0*sgm12/rij2**7+sgm6/rij2**4)
    !力 
                  dfdx(1,i)=dfdx(1,i)+factor*rijx
                  dfdx(2,i)=dfdx(2,i)+factor*rijy
                  dfdx(3,i)=dfdx(3,i)+factor*rijz
                  dfdx(1,j)=dfdx(1,j)-factor*rijx
                  dfdx(2,j)=dfdx(2,j)-factor*rijy
                  dfdx(3,j)=dfdx(3,j)-factor*rijz
    !圧力
                  p(1,1)=p(1,1)+factor*rijx*rijx
                  p(1,2)=p(1,2)+factor*rijx*rijy
                  p(1,3)=p(1,3)+factor*rijx*rijz
                  p(2,1)=p(2,1)+factor*rijy*rijx
                  p(2,2)=p(2,2)+factor*rijy*rijy
                  p(2,3)=p(2,3)+factor*rijy*rijz
                  p(3,1)=p(3,1)+factor*rijz*rijx
                  p(3,2)=p(3,2)+factor*rijz*rijy
                  p(3,3)=p(3,3)+factor*rijz*rijz
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
      write(*,*)p(1,1),p(2,2),p(3,3)
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

!逆行列と体積
      subroutine inverse_mass(h,h_inver,vol)
      implicit none
      real*8 h(3,3)
      real*8 h_inver(3,3)
      real*8 vol

      vol=h(1,1)*(h(2,2)*h(3,3)-h(2,3)*h(3,2))
     &   -h(1,2)*(h(2,1)*h(3,3)-h(2,3)*h(3,1))
     &   +h(1,3)*(h(2,1)*h(3,2)-h(2,2)*h(3,1))

      h_inver(1,1)=(h(2,2)*h(3,3)-h(2,3)*h(3,2))/vol
      h_inver(1,2)=(h(2,3)*h(3,1)-h(2,1)*h(3,3))/vol
      h_inver(1,3)=(h(2,1)*h(3,2)-h(3,1)*h(2,2))/vol
      h_inver(2,1)=(h(1,3)*h(3,2)-h(1,2)*h(3,3))/vol
      h_inver(2,2)=(h(1,1)*h(3,3)-h(1,3)*h(3,1))/vol
      h_inver(2,3)=(h(1,2)*h(3,1)-h(1,1)*h(3,2))/vol
      h_inver(3,1)=(h(1,2)*h(2,3)-h(1,3)*h(2,2))/vol
      h_inver(3,2)=(h(1,3)*h(2,1)-h(1,1)*h(2,3))/vol
      h_inver(3,3)=(h(1,1)*h(2,2)-h(1,2)*h(2,1))/vol
      end subroutine

! 圧力      
      subroutine pres(v,n,p,vol)
      implicit none
      integer n,i,j,m
      real*8 v(3,n),p(3,3)
      real*8 vol,ekin,Tk,amass,kb

      amass=40d0*1836d0
      kb=1.d0/(27.2116*11605d0)

      ekin=0d0
      do i=1,n
        ekin=ekin+0.5d0*amass*
     &   (v(1,i)**2+v(2,i)**2+v(3,i)**2)
      enddo
      Tk=ekin*2d0/(3d0*dble(n))*
     &    27.2116*11605d0
      do j=1,3
      do m=1,3
        if(j==m) then
        p(j,m)=(n*kb*Tk+p(j,m))/vol
        else
        p(j,m)=p(j,m)/vol
       endif
      enddo
      enddo
      end subroutine




