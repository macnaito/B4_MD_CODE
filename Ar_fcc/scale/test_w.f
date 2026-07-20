!　スケールした座標でのMD エネルギー保存を確認
!　lclを少し改良　parrinello_rahman+nose_hoover
!　変数のWとtauを決めたい

      program ljmd
      
      implicit real*8 (a-h,o-z)
      real*8 kb
      parameter(ev_k=11605.d0)
      parameter(hartree_ev=27.2116d0)
      parameter(kb=1.d0/ev_k/hartree_ev)
      parameter(Gp=29421.d0)
      parameter(bohr=0.5292d0)

      real*8,allocatable,dimension(:,:)::x,v,frc,s,ds
      real*8,allocatable,dimension(:)::mass
      character*8,allocatable,dimension(:)::lsp
      character*40 filename
      real*8 p(3,3),h(3,3),h_inver(3,3),sgm(3,3)
      real*8 str(3,3),e_kin(3,3)
      real*8 dh(3,3),Preg(3,3),GinvGdot(3,3)
      real*8 Treg,Q,tau,xi,gkbt,eta
      real*8 tauloop(7)
      integer c

      real*8,allocatable,dimension(:)::rec1,rec2,rec3,rec4,rec5,rec6
      real*8,allocatable,dimension(:)::gpa,smap

      maxstep=3000
      allocate(gpa(maxstep),smap(maxstep),rec1(maxstep),rec2(maxstep)
     &        ,rec3(maxstep),rec4(maxstep),rec5(maxstep),rec6(maxstep))

      tauloop=[40,50,60,70,80,90,100]

      do c=1,7
      tau=tauloop(c)*dt

      k=0
      eta=0.d0
      hmax=-1.d10
      hmin=1.d10
      pmax=-1d10
      pmin=1d10
      xi=0.d0

!　.datファイルの読み込み      
      open(10,file='fainal.dat')
       read(10,*)ntot
        allocate(x(3,ntot),v(3,ntot),frc(3,ntot),s(3,ntot),ds(3,ntot)
     &           ,mass(ntot),lsp(ntot))
       do i=1,ntot
        read(10,*)lsp(i),dummy,x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i)
        if(lsp(i)=='Ar') mass(i)=40d0*1836d0
        if(lsp(i)/='Ar') stop
       enddo
       read(10,*) h(1,1),h(2,1),h(3,1)
       read(10,*) h(1,2),h(2,2),h(3,2)
       read(10,*) h(1,3),h(2,3),h(3,3)
      close(10)
      h=h/bohr
      x(:,:)=x(:,:)/bohr    !a.u.に変換
!
      dt=41.d0*3.d0
      call vol_inverse(h,vol,h_inver,sgm)
!　スケール  
      do i=1,ntot
        s(:,i)=matmul(h_inver,x(:,i))
        ds(:,i)=matmul(h_inver,v(:,i))
      enddo 
      call pot(f,frc,s,x,ntot,h,str)

! 温度圧力の計算
      ekin=0.d0
      do i=1,ntot
       ekin=ekin+0.5d0*mass(i)*sum(v(:,i)**2)
      enddo
      tempk=ekin*2.d0/3.d0/ntot/kb
      call p_tesor(s,ds,h,dh,str,ntot,mass,vol,p,e_kin,ekin)
      gpa(1)=(p(1,1)+p(2,2)+p(3,3))/3.d0
!
      Preg=0.d0
      Preg(1,1)=1d-4/Gp; Preg(2,2)=1d-4/Gp; Preg(3,3)=1d-4/Gp
      dh=0.d0
      w=3.d9
      Treg=50.d0
 !     tau=65.d0*dt
      gkbt=3.d0*ntot*Treg*kb
      Q=gkbt*tau**2
      xi=0.d0

 !     write(*,'(A,1x,F0.2,A,4x,A,1x,F0.2,A)')
 !    &  '初期温度',tempk,'K','初期圧力',gpa(1)*Gp*1d4,'気圧'
 !     write(*,'(A,1x,F0.2,A,4x,A,1x,F0.2,A)')
 !    &  '目標温度',Treg,'K','目標圧力',Preg(1,1)*Gp*1d4,'気圧'
 !     write(*,'(A,f0.2,A,2x,A,I0,A)')
 !    & 'dt:',dt/41.d0,'fs','maxstep:',maxstep,'step'
 !     write(*,*)'温度制御nosehoover 圧力制御parirellorahman'

 !     read(*,*)

!=========================計算スタート====================
      do istep=1,maxstep

!　速度の更新
      GinvGdot=matmul(dh,h_inver)+transpose(matmul(dh,h_inver))
      do i=1,ntot
        ds(:,i)=ds(:,i)+0.5d0*dt*matmul(h_inver,frc(:,i)/mass(i))
     &                 -0.5d0*dt*matmul(GinvGdot,ds(:,i))
     &                 -0.5d0*dt*xi*ds(:,i)
      enddo
      call p_tesor(s,ds,h,dh,str,ntot,mass,vol,p,e_kin,ekin)
      dh=dh+0.5d0*dt*matmul((p-Preg),sgm)/w
      xi=xi+0.5d0*dt*(2.d0*ekin-gkbt)/Q
      eta=eta+xi*0.5d0*dt
!         
!　位置の更新+pbc
      do i=1,ntot
        s(:,i)=s(:,i)+dt*ds(:,i)
        do j=1,3
          if(s(j,i)>=1.d0) s(j,i)=s(j,i)-1.d0
          if(s(j,i)< 0.d0) s(j,i)=s(j,i)+1.d0
        enddo
      enddo
      h=h+dt*dh
! 
      x=matmul(h,s)
      call vol_inverse(h,vol,h_inver,sgm)
      call pot(f,frc,s,x,ntot,h,str)

!　速度の更新
      GinvGdot=matmul(dh,h_inver)+transpose(matmul(dh,h_inver))
      do i=1,ntot
        ds(:,i)=ds(:,i)+0.5d0*dt*matmul(h_inver,frc(:,i)/mass(i))
     &                 -0.5d0*dt*matmul(GinvGdot,ds(:,i))
     &                 -0.5d0*dt*xi*ds(:,i)
      enddo
      call p_tesor(s,ds,h,dh,str,ntot,mass,vol,p,e_kin,ekin)
      dh=dh+0.5d0*dt*matmul((p-Preg),sgm)/w
      xi=xi+0.5d0*dt*(2.d0*ekin-gkbt)/Q
      eta=eta+xi*0.5d0*dt
!
!　温度・圧力の計算
      tempk=ekin*2.d0/3.d0/ntot/kb
      gpa(istep)=(p(1,1)+p(2,2)+p(3,3))/3.d0
      if(istep>400) then
       pmax=max(pmax,gpa(istep));pmin=min(pmin,gpa(istep))
      endif


!　xyzファイルの出力
      filename='out000.xyz'
      if(mod(istep,maxstep/100).eq.0)then
      k=k+1
      do i=1,ntot
        x(:,i)=matmul(h,s(:,i))
      enddo
      write(filename(4:6),'(i3.3)')k
      open(11,file=filename)
      write(11,*)ntot
      write(11,'(a,3(3e15.7),a,a)')
     &       'Lattice="',h(1,1)*bohr,h(2,1)*bohr,h(3,1)*bohr,
     &                   h(1,2)*bohr,h(2,2)*bohr,h(3,2)*bohr,
     &                   h(1,3)*bohr,h(2,3)*bohr,h(3,3)*bohr,'" ',
     &       'Properties=species:S:1:id:I:1:pos:R:3:tempK:R:1'
      do m=1,ntot
      write(11,'(a2,1x,i5,4e15.7)')lsp(m),m,
     &       x(1,m)*bohr,x(2,m)*bohr,x(3,m)*bohr,0d0
      enddo
      close(11)
      endif
!
      hamil=ekin+f+0.5d0*w*sum(dh**2)
     &    +vol*(Preg(1,1)+Preg(2,2)+Preg(3,3))/3.d0
     &    +0.5d0*Q*xi**2+gkbt*eta
      hmax=max(hmax,hamil);hmin=min(hmin,hamil)
 !     write(*,*)istep,tempk,gpa(istep)*Gp,hamil
      
      rec1(istep)=tempk
      rec2(istep)=f
      rec3(istep)=ekin
      rec4(istep)=0.5d0*w*sum(dh**2)
      rec5(istep)=vol*Preg(1,1)
      rec6(istep)=hamil
      
      enddo
!=====================loopここまで==============================

      write(*,*)c
      write(*,"(5x,F0.4,A)")
     &     abs((hmax-hmin)/hamil*100),'%'
      write(*,*)pmax*Gp,pmin*Gp,(pmax-pmin)*Gp
      deallocate(x,v,frc,s,ds,mass,lsp)

      enddo

      

!　グラフ
      smap=0.d0
      call sma(gpa,smap,maxstep)
      open(12,file='testw.dat')
      do i=1,maxstep
        write(12,*)rec1(i),rec2(i),rec3(i),rec4(i),rec5(i),rec6(i)
     &             ,gpa(i)*gp,smap(i)*Gp
      enddo
      close(12)  
!

      end program


!　体積と逆行列
      subroutine vol_inverse(h,vol,h_inver,sgm)
      implicit real*8 (a-h,o-z)
      real*8 h(3,3),h_inver(3,3),sgm(3,3)
       
      vol=h(1,1)*h(2,2)*h(3,3)+h(1,2)*h(2,3)*h(3,1)
     &   +h(1,3)*h(2,1)*h(3,2)-h(1,3)*h(2,2)*h(3,1)
     &   -h(1,2)*h(2,1)*h(3,3)-h(1,1)*h(2,3)*h(3,2)  
 
      sgm(1,1)=(h(2,2)*h(3,3)-h(2,3)*h(3,2))
      sgm(1,2)=(h(1,3)*h(3,2)-h(1,2)*h(3,3))
      sgm(1,3)=(h(1,2)*h(2,3)-h(1,3)*h(2,2))
      sgm(2,1)=(h(2,3)*h(3,1)-h(2,1)*h(3,3))
      sgm(2,2)=(h(1,1)*h(3,3)-h(1,3)*h(3,1))
      sgm(2,3)=(h(1,3)*h(2,1)-h(1,1)*h(2,3))
      sgm(3,1)=(h(2,1)*h(3,2)-h(2,2)*h(3,1))
      sgm(3,2)=(h(1,2)*h(3,1)-h(1,1)*h(3,2))
      sgm(3,3)=(h(1,1)*h(2,2)-h(1,2)*h(2,1))

      h_inver=sgm/vol

      return
      end
!
!　call pot
      subroutine pot(f,frc,s,x,ntot,h,str)
      implicit real*8 (a-h,o-z)
      integer hx_lc,hy_lc,hz_lc,hxyz_lc,hyz_lc,pea
      real*8 frc(3,ntot),x(3,ntot),h(3,3),str(3,3),rij(3)
      real*8 hx,hy,hz,t(3,8),s(3,ntot)
      parameter(sgm=3.4d0/0.5292d0)  
      parameter(eps=120d0/11605d0/27.2116d0)
      parameter(sgm12=sgm**12,sgm6=sgm**6)
      parameter(cutoff=2.5d0*sgm)
      integer,allocatable,dimension(:)::lshd,lscl

      frc=0.d0
      str=0.d0
      pea=0

 !シミュレーションボックスが入る立方体の各頂点と最小最大座標、立方体の各辺の長さ
      t(:,1)=(/0.d0,0.d0,0.d0/)
      t(:,2)=h(:,1)
      t(:,3)=h(:,1)+h(:,2)
      t(:,4)=h(:,2)
      t(:,5)=h(:,3)
      t(:,6)=h(:,1)+h(:,3)
      t(:,7)=h(:,1)+h(:,2)+h(:,3)
      t(:,8)=h(:,2)+h(:,3)
      xmin=minval(t(1,:))
      xmax=maxval(t(1,:))
      ymin=minval(t(2,:))
      ymax=maxval(t(2,:))
      zmin=minval(t(3,:))
      zmax=maxval(t(3,:))
      hx=xmax-xmin
      hy=ymax-ymin
      hz=zmax-zmin

      hx_lc=int(hx/cutoff) !x座標のセル数
      hy_lc=int(hy/cutoff)
      hz_lc=int(hz/cutoff)
      hyz_lc=hy_lc*hz_lc     !セルのyz平面の数
      hxyz_lc=hyz_lc*hx_lc   !総セル数
      hx_cell=hx/hx_lc     !各方向のセルの長さ
      hy_cell=hy/hy_lc
      hz_cell=hz/hz_lc
      
      allocate(lshd(hxyz_lc),lscl(ntot))
      lshd=0
      do i=1,ntot
        rx=x(1,i)-xmin
        ry=x(2,i)-ymin
        rz=x(3,i)-zmin
        mx=int(rx/hx_cell)
        my=int(ry/hy_cell)
        mz=int(rz/hz_cell)
        m=mx*hyz_lc+my*hz_lc+mz+1
        lscl(i)=lshd(m)
        lshd(m)=i
      enddo

      kuxmax=ceiling(cutoff/hx_cell)
      kuymax=ceiling(cutoff/hy_cell)
      kuzmax=ceiling(cutoff/hz_cell)
      f=0d0
      
      do mz=0,hz_lc-1
      do my=0,hy_lc-1
      do mx=0,hx_lc-1
        m=mx*hyz_lc+my*hz_lc+mz+1 !3つのdoloopで全てのセルを走査
        if(lshd(m)==0) cycle !空のセルはスキップ
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
                rij=rij-dnint(rij)   
                rij=matmul(h,rij)        
                rij2=rij(1)**2+rij(2)**2+rij(3)**2
                if (rij2<(cutoff)**2) then
                  pea=pea+1
                  f=f+4d0*eps*(sgm12/rij2**6-sgm6/rij2**3)
                  factor=4d0*eps*
     &            (-12d0*sgm12/rij2**7+6d0*sgm6/rij2**4)
                  frc(:,i)=frc(:,i)-factor*rij(:)
                  frc(:,j)=frc(:,j)+factor*rij(:)
                  do k=1,3
                  do l=1,3
                   str(k,l)=str(k,l)-factor*rij(k)*rij(l)
                  enddo
                  enddo
                endif
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
      return
      end
!
! subrutine ishift  
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
!
!　単純移動平均 simple-moving-average
      subroutine sma(gpa,smap,maxstep)
      implicit none
      real*8 gpa(maxstep),smap(maxstep),sum
      integer maxstep,n,i
      parameter(n=500)

      sum=0.d0
      do i=1,maxstep
        if(i<=n) then
        sum=sum+gpa(i)
        elseif(i>n) then
        sum=sum+gpa(i)-gpa(i-n)
        smap(i)=sum/n
        endif
      enddo

      return
      end
!
!　圧力テンソル
      subroutine p_tesor(s,ds,h,dh,str,ntot,mass,vol,p,e_kin,ekin)
      implicit none
      integer ntot,i,k,l
      real*8 s(3,ntot),ds(3,ntot),str(3,3),p(3,3),h(3,3)
      real*8 e_kin(3,3),mass(ntot),vol,vi(3),dh(3,3),ekin
            
      e_kin=0.d0
      ekin=0.d0
      do i=1,ntot
        vi=matmul(h,ds(:,i))+matmul(dh,s(:,i))
       do k=1,3
       do l=1,3  
       e_kin(k,l)=e_kin(k,l)+mass(i)*vi(k)*vi(l)
       enddo
       enddo
      enddo
      ekin=(e_kin(1,1)+e_kin(2,2)+e_kin(3,3))*0.5d0
      p=(e_kin+str)/vol
      return
      end
!
!　記録
      subroutine kiroku(h,ntot,bohr,x,v)
      implicit none
      real*8 h(3,3),x(3,ntot),v(3,ntot),bohr
      integer i,ntot
      
      open(13,file='fainal.dat')
      write(13,*)ntot
      do i=1,ntot
        write(13,'(a,1x,i5,6e15.7)') 'Ar',i,
     &       x(1,i)*bohr,x(2,i)*bohr,x(3,i)*bohr,
     &       v(1,i),v(2,i),v(3,i)
      enddo
      write(13,'(3e24.15)')h(1,1)*bohr,h(2,1)*bohr,h(3,1)*bohr
      write(13,'(3e24.15)')h(1,2)*bohr,h(2,2)*bohr,h(3,2)*bohr
      write(13,'(3e24.15)')h(1,3)*bohr,h(2,3)*bohr,h(3,3)*bohr
      close(13)
      end
!