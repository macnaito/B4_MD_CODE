      program pbc_lv_pr
! PBC lcl lengevine（BAOAB) Parrinello-Rahman
      implicit none
      integer nmax
      parameter(nmax=100000)
      integer i,j,k,m,n,maxstep,itemp,istep
      real*8 x(3,nmax),v(3,nmax),dfdx(3,nmax),f,dt
      real*8 h(3,3),vol,amass,bohr
      character*2 lsp(nmax)
      character*40 filename2
      parameter(amass=40d0*1836d0)
      parameter(bohr=0.5292d0)
      real*8 W,sgm(3,3),tk,kin,rmin
      real*8 virial(3,3),p(3,3),kb,pext(3,3)
      real*8 s(3,nmax),vh(3,3),h_inver(3,3)
      real*8 vs(3,nmax),mt(3,3),mt_inver(3,3)
      real*8 Treg,tau,gkbt,Q,xi,G(3,3)
      real*8 smap
 
      
! time_step      
      parameter(maxstep=10000)
      real*8 gpa(maxstep)
      real*8 t(maxstep)
      real*8 hx(maxstep)
      real*8 smagpa(maxstep)
    

      dt=41d0*5d0
      kb=1.d0/(27.2116*11605d0)
      W=5.d5
      xi=0d0

      pext=0d0
      pext(1,1)=1.01325d-4/29421d0
      pext(2,2)=1.01325d-4/29421d0
      pext(3,3)=1.01325d-4/29421d0
      filename2='out000.xyz'
      k=0
      vh=0d0

!ファイルの読みこみ      
      open(10,file='lv50_bs1d-4.dat')
      read(10,*)n
      do i=1,n
        read(10,*)lsp(i),itemp,
     &  x(1,i),x(2,i),x(3,i),
     &  v(1,i),v(2,i),v(3,i)

        x(1,i)=x(1,i)/bohr
        x(2,i)=x(2,i)/bohr
        x(3,i)=x(3,i)/bohr
      enddo
      read(10,*)h(1,1),h(2,1),h(3,1)
      read(10,*)h(1,2),h(2,2),h(3,2)
      read(10,*)h(1,3),h(2,3),h(3,3)
      h=h/bohr
      close(10)
!ファイルの読み込み

      Treg=20d0
      tau=10d0*dt
      gkbt=(3.d0*dble(n))*Treg*kb
      Q=gkbt*tau**2


!計算start    

!温度
      kin=0d0
      do j=1,n
          kin=kin+0.5d0*amass*
     &     (v(1,j)**2+v(2,j)**2+v(3,j)**2)
      enddo
      tk=kin*2d0/(3d0*dble(n))*
     &     27.2116*11605d0
!温度
    
      call inverse_mass(h,h_inver,vol,mt,mt_inver,sgm)
      do i=1,n
       s(:,i)=matmul(h_inver,x(:,i)) 
       vs(:,i)=matmul(h_inver,v(:,i))
      enddo
  
      
  
      call pot (f,dfdx,s,n,h,virial,rmin) ! mt:metoric tensor
      call inverse_mass(h,h_inver,vol,mt,mt_inver,sgm)
      call press(v,n,vol,p,virial)
      write(*,*)p(1,1)*29421d0
      

! ここからloop
      do istep=1,maxstep

      Treg=50d0
      if (istep<=5000) then
         Treg=Treg+0.002d0*istep
        else
         Treg=Treg+10d0
        endif
        gkbt=(3.d0*dble(n))*Treg*kb
        Q=gkbt*tau**2

        kin=0d0
        do j=1,n
          kin=kin+0.5d0*amass*
     &     (v(1,j)**2+v(2,j)**2+v(3,j)**2)
        enddo
        xi=xi+0.5d0*dt*(2.d0*kin-gkbt)/Q


        v=matmul(h,vs)
        G=matmul(vh,h_inver)
        do i=1,n
          v(:,i)=v(:,i)+0.5d0*dt*(dfdx(:,i)/amass)
     &             -0.5d0*dt*xi*v(:,i)    !nose-hoover
     &             -0.5d0*dt*matmul(G,v(:,i))  !
          vs(:,i)=matmul(h_inver,v(:,i))
        enddo
        vh=vh+0.5d0*dt*(matmul((p-pext),sgm)/W)

  
!位置の更新        
        do i=1,n
         s(:,i)=s(:,i)+dt*vs(:,i)
        enddo
        h=h+dt*vh

!pbc
        do j=1,3
          do i=1,n
           if(s(j,i) >= 1.d0) s(j,i)=s(j,i)-1.d0
           if(s(j,i) < 0.d0)  s(j,i)=s(j,i)+1.d0
          enddo
        enddo
!pbc


        x=matmul(h,s)
        v=matmul(h,vs)

        call pot(f,dfdx,s,n,h,virial,rmin)
        call inverse_mass(h,h_inver,vol,mt,mt_inver,sgm)
        call press(v,n,vol,p,virial)

        kin=0d0
        do j=1,n
          kin=kin+0.5d0*amass*
     &     (v(1,j)**2+v(2,j)**2+v(3,j)**2)
          tk=kin*2d0/(3d0*dble(n))*
     &     27.2116*11605d0
        enddo
        xi=xi+0.5d0*dt*(2.d0*kin-gkbt)/Q


        v=matmul(h,vs)
        do i=1,n
         v(:,i)=v(:,i)+0.5d0*dt*(dfdx(:,i)/amass)
     &                -0.5d0*dt*xi*v(:,i)
     &                -0.5d0*dt*matmul(G,v(:,i))
         vs(:,i)=matmul(h_inver,v(:,i))
        enddo
        vh=vh+0.5d0*dt*matmul((p-pext),sgm)/W    
     
!ファイルの書き出し
        if(mod(istep,100).eq.0)then
          k=k+1
          write(filename2(4:6),'(i3.3)')k

          open(11,file=filename2)
          write(11,*)n
          write(11,'(a,3(3e15.7),a,a)')
     &       'Lattice="',h(1,1)*bohr,h(2,1)*bohr,h(3,1)*bohr,
     &                   h(1,2)*bohr,h(2,2)*bohr,h(3,2)*bohr,
     &                   h(1,3)*bohr,h(2,3)*bohr,h(3,3)*bohr,'" ',
     &       'Properties=species:S:1:id:I:1:pos:R:3:tempK:R:1'
           do m=1,n
            write(11,'(a2,1x,i5,4e15.7)')lsp(m),m,
     &       x(1,m)*bohr,x(2,m)*bohr,x(3,m)*bohr,0d0
          enddo
          close(11)
        endif
!ファイルへの書き出し 

 
        write(*,*) istep,p(1,1)*29421d0,h(1,1),tk
        hx(istep)=(h(1,1)+h(2,2)+h(3,3))/3.d0
        t(istep)=tk
        gpa(istep)=(p(1,1)+p(2,2)+p(3,3))*29421d0

!        call simple_moving_average(istep,gpa,smap)
!        write(*,*)smap
!        smagpa(istep)=smap


      enddo

!グラフ      
      open (12,file='t.dat')
      do i=1,maxstep
        write(12,*) gpa(i),t(i),hx(i)
      enddo
      close(12)
!グラフ

!kiroku
      call kiroku(h,n,bohr,x,v)
!kiroku

      end program pbc_lv_pr


!linked cell
      subroutine pot(f,dfdx,s,n,h,virial,rmin)
      implicit none
      integer mx,my,mz,hx_lc,hy_lc,hz_lc,m1,m
      real*8 dfdx(3,n),f,factor,h(3,3),virial(3,3)
      real*8 s(3,n)
      real*8 sgm,eps,sgm12,sgm6,cutoff,rij2
      integer hyz_lc,hxyz_lc,i,m1x,m1y,m1z
      integer ishiftx,ishifty,ishiftz,j,k,l,n
      integer kuxmax,kuymax,kuzmax,kux,kuy,kuz
      parameter(sgm=3.4d0/0.5292d0)
      parameter(eps=120d0/11605d0/27.2116d0)
      parameter(sgm12=sgm**12,sgm6=sgm**6)
      parameter(cutoff=2.5d0*sgm)
      integer,allocatable,dimension(:):: lshd,lscl
      real*8 rij(3),rmin,hx_cell,hy_cell,hz_cell
      integer counter,counterout

      counter=0
      counterout=0

      dfdx=0d0
      rmin=1000d0
      virial=0d0

       hx_lc=ceiling(sqrt(dot_product(h(:,1),h(:,1)))/cutoff)-1 !x座標のセル数
       hy_lc=ceiling(sqrt(dot_product(h(:,2),h(:,2)))/cutoff)-1
       hz_lc=ceiling(sqrt(dot_product(h(:,3),h(:,3)))/cutoff)-1
       hyz_lc=hy_lc*hz_lc     !セルのyz平面の数
       hxyz_lc=hyz_lc*hx_lc     !セルの総数
       hx_cell=1.d0/hx_lc     !各方向のセルの長さ
       hy_cell=1.d0/hy_lc
       hz_cell=1.d0/hz_lc
       allocate(lshd(hxyz_lc),lscl(n))
       lshd=0
       do i=1,n
         mx=int(s(1,i)/hx_cell)
         my=int(s(2,i)/hy_cell)
         mz=int(s(3,i)/hz_cell)
         mx=min(max(mx,0),hx_lc-1)
         my=min(max(my,0),hy_lc-1)
         mz=min(max(mz,0),hz_lc-1)
         m=mx*hyz_lc+my*hz_lc+mz+1
         lscl(i)=lshd(m)
         lshd(m)=i
       enddo

       kuxmax=1
       kuymax=1
       kuzmax=1
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
                rij(:)=s(:,i)-s(:,j)
                rij(:)=rij(:)-dnint(rij(:))
                rij=matmul(h,rij)
                rij2=rij(1)**2+rij(2)**2+rij(3)**2
                if (rij2<cutoff**2) then
 !                 counter=counter+1
                  f=f+4d0*eps*(sgm12/rij2**6-sgm6/rij2**3)
                  factor=4d0*eps*
     &            (-12d0*sgm12/rij2**7+6d0*sgm6/rij2**4)
                  dfdx(1,i)=dfdx(1,i)-factor*rij(1)
                  dfdx(2,i)=dfdx(2,i)-factor*rij(2)
                  dfdx(3,i)=dfdx(3,i)-factor*rij(3)
                  dfdx(1,j)=dfdx(1,j)+factor*rij(1)
                  dfdx(2,j)=dfdx(2,j)+factor*rij(2)
                  dfdx(3,j)=dfdx(3,j)+factor*rij(3)

                  do k=1,3
                  do l=1,3
                   virial(k,l)=virial(k,l)-factor*rij(k)*rij(l)
                  enddo
                  enddo
                else
 !                  counterout=counterout+1
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
 !     write(*,*)counter,counterout,dble(counter)/(counter+counterout)
      return
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
      subroutine inverse_mass(h,h_inver,vol,mt,mt_inver,sgm)
      implicit none
      real*8 h(3,3)
      real*8 h_inver(3,3)
      real*8 mt(3,3),mt_inver(3,3),mt_vol
      real*8 vol,sgm(3,3)

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

      mt=matmul(transpose(h),h)
      mt_vol=mt(1,1)*(mt(2,2)*mt(3,3)-mt(2,3)*mt(3,2))
     &   -mt(1,2)*(mt(2,1)*mt(3,3)-mt(2,3)*mt(3,1))
     &   +mt(1,3)*(mt(2,1)*mt(3,2)-mt(2,2)*mt(3,1))
      mt_inver(1,1)=(mt(2,2)*mt(3,3)-mt(2,3)*mt(3,2))/mt_vol
      mt_inver(1,2)=(mt(2,3)*mt(3,1)-mt(2,1)*mt(3,3))/mt_vol
      mt_inver(1,3)=(mt(2,1)*mt(3,2)-mt(3,1)*mt(2,2))/mt_vol
      mt_inver(2,1)=(mt(1,3)*mt(3,2)-mt(1,2)*mt(3,3))/mt_vol
      mt_inver(2,2)=(mt(1,1)*mt(3,3)-mt(1,3)*mt(3,1))/mt_vol
      mt_inver(2,3)=(mt(1,2)*mt(3,1)-mt(1,1)*mt(3,2))/mt_vol
      mt_inver(3,1)=(mt(1,2)*mt(2,3)-mt(1,3)*mt(2,2))/mt_vol
      mt_inver(3,2)=(mt(1,3)*mt(2,1)-mt(1,1)*mt(2,3))/mt_vol
      mt_inver(3,3)=(mt(1,1)*mt(2,2)-mt(1,2)*mt(2,1))/mt_vol

      sgm(1,1)=(h(2,2)*h(3,3)-h(2,3)*h(3,2))
      sgm(1,2)=(h(1,3)*h(3,2)-h(1,2)*h(3,3))
      sgm(1,3)=(h(1,2)*h(2,3)-h(1,3)*h(2,2))
      sgm(2,1)=(h(2,3)*h(3,1)-h(2,1)*h(3,3))
      sgm(2,2)=(h(1,1)*h(3,3)-h(1,3)*h(3,1))
      sgm(2,3)=(h(1,3)*h(2,1)-h(1,1)*h(2,3))
      sgm(3,1)=(h(2,1)*h(3,2)-h(2,2)*h(3,1))
      sgm(3,2)=(h(1,2)*h(3,1)-h(1,1)*h(3,2))
      sgm(3,3)=(h(1,1)*h(2,2)-h(1,2)*h(2,1))

      return
      end subroutine

! call press      
      subroutine press(v,n,vol,p,virial)
      implicit none
      integer n,l,m
      real*8 v(3,n),p(3,3)
      real*8 vol
      real*8 amass,a(3,3),virial(3,3)

      parameter(amass=40d0*1836d0)

      a=amass*matmul(v,transpose(v))
      do l=1,3
      do m=1,3
       p(l,m)=(a(l,m)+virial(l,m))/vol
      enddo
      enddo

      return
      end subroutine

! kiroku
      subroutine kiroku(h,n,bohr,x,v)
      implicit none
      real*8 h(3,3),bohr,x(3,n),v(3,n)
      integer i,n

      open(10,file='nh50-60_pr1d-4.dat')
      write(10,*)n
      do i=1,n
       write(10,'(a,1x,i5,6e15.7)') 'Ar',i,
     &       x(1,i)*bohr,x(2,i)*bohr,x(3,i)*bohr,
     &       v(1,i),v(2,i),v(3,i)
      enddo
      write(10,'(3e24.15)')h(1,1)*bohr,h(2,1)*bohr,h(3,1)*bohr
      write(10,'(3e24.15)')h(1,2)*bohr,h(2,2)*bohr,h(3,2)*bohr
      write(10,'(3e24.15)')h(1,3)*bohr,h(2,3)*bohr,h(3,3)*bohr
      close(10)

      return
      end 

!      subroutine simple_moving_average(istep,gpa,smap)
!      implicit none
!      real*8 gpa(istep),smap
!      integer istep

!       if (istep<=)
       
!
!      return
!      end