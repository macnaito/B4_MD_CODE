!　スケールした座標でのMD 圧力制御なし　エネルギー保存を確認
!　速度こうしんはvで行う　vcl　温度制御nose-hoover
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
      real*8 p(3,3),h(3,3),h_inver(3,3)
      real*8 G(3,3),e_kin(3,3),str(3,3)
      integer,allocatable,dimension(:)::nei_counter
      integer,allocatable,dimension(:,:)::nei_list

      real*8,allocatable,dimension(:)::rec1,rec2,rec3,rec4,rec5,rec6
      real*8,allocatable,dimension(:)::gpa,smap
      real*8 Treg,Q,tau,xi,gkbt,eta,hamil

      maxstep=6000
      allocate(gpa(maxstep),smap(maxstep),rec1(maxstep),rec2(maxstep)
     &        ,rec3(maxstep),rec4(maxstep),rec5(maxstep),rec6(maxstep))

       k=0
       dt=41.d0*1.d0
       eta=0.d0
       hmax=-1.d10
       hmin=1.d10

!　.datファイルの読み込み      
      open(10,file='fainal.dat')
!
       read(10,*)ntot
        allocate(x(3,ntot),v(3,ntot),frc(3,ntot),s(3,ntot),ds(3,ntot),
     &              mass(ntot),lsp(ntot))
       do i=1,ntot
        read(10,*)lsp(i),dammy,x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i)
        if(lsp(i)=='Ar') mass(i)=40d0*1836d0
        if(lsp(i)/='Ar') stop
       enddo
       read(10,*) h(1,1),h(2,1),h(3,1)
       read(10,*) h(1,2),h(2,2),h(3,2)
       read(10,*) h(1,3),h(2,3),h(3,3)
      close(10)
      h=h/bohr
      x(:,:)=x(:,:)/bohr    !a.u.に変換
      allocate(nei_counter(ntot),nei_list(200,ntot))
!

      call vol_inverse(h,vol,h_inver)

!　スケール  
      do i=1,ntot
        s(:,i)=matmul(h_inver,x(:,i))
        ds(:,i)=matmul(h_inver,v(:,i))
      enddo 

! 温度の計算
      call cpu_time(t1)
      ekin=0.d0
      do i=1,ntot
       ekin=ekin+0.5d0*mass(i)*sum(v(:,i)**2)
      enddo
      tempk=ekin*2.d0/3.d0/ntot/kb
      call cpu_time(t2)
      write(*,*)tempk,t2-t1

      call cpu_time(t1)
      ekin=0.d0
      do i=1,ntot
       ekin=ekin+0.5d0*mass(i)
     &  *(v(1,i)**2+v(2,i)**2+v(3,i)**2)
      enddo
      tempk=ekin*2.d0/3.d0/ntot/kb
      call cpu_time(t2)
      write(*,*)tempk,t2-t1
 
      call cpu_time(t1)
      ekin=0.d0
      G=matmul(transpose(h),h)
      do i=1,ntot
       ekin=ekin+0.5d0*mass(i)*dot_product(ds(:,i),matmul(G,ds(:,i)))
      enddo
      tempk=ekin*2.d0/3.d0/ntot/kb
      call cpu_time(t2)
      write(*,*)tempk,t2-t1

      call cpu_time(t1)
      e_kin=0.d0
       do i=1,ntot
        do l=1,3
        do m=1,3
          e_kin(l,m)=e_kin(l,m)+0.5d0*mass(i)*v(l,i)*v(m,i)
        enddo
        enddo 
      enddo
      tempk=(e_kin(1,1)+e_kin(2,2)+e_kin(3,3))*2.d0/3.d0/ntot/kb
      call cpu_time(t2)
      write(*,*)tempk,t2-t1

      call cpu_time(t1)
      e_kin=0.d0
      G=matmul(transpose(h),h)
       do i=1,ntot
        do l=1,3
        do m=1,3
          e_kin(l,m)=e_kin(l,m)+0.5d0*mass(i)
     &              *ds(l,i)*dot_product(G(m,:),ds(:,i))
        enddo
        enddo 
      enddo
      tempk=(e_kin(1,1)+e_kin(2,2)+e_kin(3,3))*2.d0/3.d0/ntot/kb
      call cpu_time(t2)
      write(*,*)tempk,t2-t1
!
      call vcl(h,s,ntot,nei_counter,nei_list)
      call pot(f,frc,s,ntot,h,str,nei_counter,nei_list)

      Treg=50.d0
      tau=75.d0*dt
      gkbt=(3.d0*ntot*Treg*kb)
      Q=gkbt*tau**2
      xi=0.5d0*dt*(2.d0*ekin-gkbt)/Q

!=========================計算スタート====================
      do istep=1,maxstep

!　速度の更新
      ekin=0.d0
      do i=1,ntot
        v(:,i)=v(:,i)+0.5d0*dt*frc(:,i)/mass(i)
     &               -0.5d0*dt*xi*v(:,i)
        ekin=ekin+0.5d0*mass(i)*sum(v(:,i)**2)
        ds(:,i)=matmul(h_inver,v(:,i))
      enddo
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
! 
!     call vcl(h,s,ntot,nei_counter,nei_list)
      call pot(f,frc,s,ntot,h,str,nei_counter,nei_list)

!　速度の更新
      ekin=0.d0  
      do i=1,ntot
        v(:,i)=v(:,i)+0.5d0*dt*frc(:,i)/mass(i)
     &               -0.5d0*dt*xi*v(:,i)
        ekin=ekin+0.5d0*mass(i)*sum(v(:,i)**2)
        ds(:,i)=matmul(h_inver,v(:,i))
      enddo
      xi=xi+0.5d0*dt*(2.d0*ekin-gkbt)/Q
      eta=eta+xi*0.5d0*dt
! 

      tempk=ekin*2.d0/3.d0/ntot/kb
      call pressure(istep,maxstep,v,ntot,str,vol,p,gpa,smap)

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

  !    write(*,*)istep,tempk,gpa(istep),smap(istep),f+ekin
      hamil=ekin+f+0.5d0*Q*xi**2+gkbt*eta
      hmax=max(hmax,hamil)
      hmin=min(hmin,hamil)
      write(*,*)istep,tempk,ekin+f,gpa(istep),hamil
      
      rec1(istep)=tempk
      rec2(istep)=f
      rec3(istep)=ekin
      rec4(istep)=hamil
      rec5(istep)=nan
      rec6(istep)=nan
      
      enddo

!=====================loopここまで==============================

!　グラフ
      open(12,file='t.dat')
      do i=1,maxstep
        write(12,*)rec1(i),rec2(i),rec3(i),rec4(i),rec5(i),rec6(i)
     &             ,gpa(i),smap(i)
      enddo
      close(12)  
!
!　記録
      call kiroku(h,ntot,bohr,x,v)
!
      write(*,*)hmax,hmin,abs((hmax-hmin)/hamil*100)

      end program





!　体積と逆行列
      subroutine vol_inverse(h,vol,h_inver)
      implicit real*8 (a-h,o-z)
      real*8 h(3,3),h_inver(3,3)
       
      vol=h(1,1)*h(2,2)*h(3,3)+h(1,2)*h(2,3)*h(3,1)
     &   +h(1,3)*h(2,1)*h(3,2)-h(1,3)*h(2,2)*h(3,1)
     &   -h(1,2)*h(2,1)*h(3,3)-h(1,1)*h(2,3)*h(3,2)  
 
      h_inver(1,1)=(h(2,2)*h(3,3)-h(2,3)*h(3,2))/vol
      h_inver(1,2)=(h(1,3)*h(3,2)-h(1,2)*h(3,3))/vol
      h_inver(1,3)=(h(1,2)*h(2,3)-h(1,3)*h(2,2))/vol
      h_inver(2,1)=(h(2,3)*h(3,1)-h(2,1)*h(3,3))/vol
      h_inver(2,2)=(h(1,1)*h(3,3)-h(1,3)*h(3,1))/vol
      h_inver(2,3)=(h(1,3)*h(2,1)-h(1,1)*h(2,3))/vol
      h_inver(3,1)=(h(2,1)*h(3,2)-h(2,2)*h(3,1))/vol
      h_inver(3,2)=(h(1,2)*h(3,1)-h(1,1)*h(3,2))/vol
      h_inver(3,3)=(h(1,1)*h(2,2)-h(1,2)*h(2,1))/vol

      return
      end
!
!　listの作成
      subroutine vcl(h,s,ntot,nei_counter,nei_list)
      implicit real*8 (a-h,o-z)
      integer hx_lc,hy_lc,hz_lc,hxyz_lc,hyz_lc
      real*8 h(3,3),s(3,ntot),rij(3)
      parameter(sgm=3.4d0/0.5292d0)
      parameter(cutoff=2.5d0*sgm)
      integer,allocatable,dimension(:)::lshd,lscl
      integer nei_counter(ntot),nei_list(200,ntot)


      nei_counter=0
      hx_lc=ceiling(sqrt(dot_product(h(:,1),h(:,1)))/cutoff)-1 !x座標のセル数
      hy_lc=ceiling(sqrt(dot_product(h(:,2),h(:,2)))/cutoff)-1
      hz_lc=ceiling(sqrt(dot_product(h(:,3),h(:,3)))/cutoff)-1
      hyz_lc=hy_lc*hz_lc     !セルのyz平面の数
      hxyz_lc=hyz_lc*hx_lc
      hx_cell=1.d0/hx_lc     !各方向のセルの長さ
      hy_cell=1.d0/hy_lc
      hz_cell=1.d0/hz_lc
      allocate(lshd(hxyz_lc),lscl(ntot))
      lshd=0
      do i=1,ntot
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
                rij(:)=rij(:)-anint(rij(:))
                rij=matmul(h,rij)
                rij2=rij(1)**2+rij(2)**2+rij(3)**2
                if (rij2<(cutoff*1.2d0)**2) then
                 nei_counter(i)=nei_counter(i)+1
                 nei_counter(j)=nei_counter(j)+1
                 nei_list(nei_counter(i),i)=j
                 nei_list(nei_counter(j),j)=i 
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

! call pot
      subroutine pot(f,frc,s,ntot,h,str,nei_counter,nei_list)
      implicit real*8 (a-h,o-z)
      real*8 frc(3,ntot),s(3,ntot),h(3,3),str(3,3),rij(3)
      parameter(sgm=3.4d0/0.5292d0)  
      parameter(eps=120d0/11605d0/27.2116d0)
      parameter(sgm12=sgm**12,sgm6=sgm**6)
      parameter(cutoff=2.5d0*sgm)
      integer nei_counter(ntot),nei_list(200,ntot),a,b

      frc=0.d0
      str=0.d0
      f=0.d0

      do a=1,ntot
      do b=1,nei_counter(a)
       j=nei_list(b,a)
       if (a<j)then
       rij=s(:,a)-s(:,j)
       rij=rij-anint(rij)
       rij=matmul(h,rij)
       rij2=dot_product(rij,rij)
       if (rij2<cutoff**2)then
        f=f+4d0*eps*(sgm12/rij2**6-sgm6/rij2**3)
        factor=4d0*eps*
     &       (-12d0*sgm12/rij2**7+6d0*sgm6/rij2**4)
        frc(:,a)=frc(:,a)-factor*rij(:)
        frc(:,j)=frc(:,j)+factor*rij(:)
        do k=1,3
        do l=1,3
         str(k,l)=str(k,l)-factor*rij(k)*rij(l)
        enddo
        enddo
       endif
       endif
      enddo
      enddo

      return
      end


!
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
!
!　瞬間の圧力と平均の圧力
      subroutine pressure(istep,maxstep,v,ntot,str,vol,p,gpa,smap)
      implicit real*8 (a-g,o-z)
      real*8 v(3,ntot),str(3,3),p(3,3),e_kin(3,3)
      real*8 smap(maxstep),gpa(maxstep)
      parameter(mass=40d0*1836d0)
      parameter(n=500)
      save sum
      data sum/0.d0/

      e_kin=mass*matmul(v,transpose(v))
      do l=1,3
      do m=1,3
        p(l,m)=(e_kin(l,m)+str(l,m))/vol
      enddo
      enddo 
      gpa(istep)=(p(1,1)+p(2,2)+p(3,3))/3.d0*29421.d0

      if(istep<=n) then
        sum=sum+gpa(istep)
        smap(istep)=sum/istep
      else 
        sum=sum-gpa(istep-n)+gpa(istep)
        smap(istep)=sum/n
      endif

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

