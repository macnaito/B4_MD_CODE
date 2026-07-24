!　スケールした座標でのMD エネルギー保存を確認
!　lcl parrinello-rahman
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
      real*8 dh(3,3),Preg(3,3),GinvGdot(3,3),GinvGdot_(3,3),G(3,3)

      real*8,allocatable,dimension(:)::rec1,rec2,rec3,rec4,rec5,rec6
      real*8,allocatable,dimension(:)::gpa,smap

      maxstep=2000
      allocate(gpa(maxstep),smap(maxstep),rec1(maxstep),rec2(maxstep)
     &        ,rec3(maxstep),rec4(maxstep),rec5(maxstep),rec6(maxstep))

       k=0

!　.datファイルの読み込み      
      open(10,file='fainal.dat')
       read(10,*)ntot
        allocate(x(3,ntot),v(3,ntot),frc(3,ntot),s(3,ntot),ds(3,ntot)
     &           ,mass(ntot),lsp(ntot))
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
!
      dt=41*1.d0

      call vol_inverse(h,vol,h_inver,sgm)

!　スケール  
      do i=1,ntot
        s(:,i)=matmul(h_inver,x(:,i))
        ds(:,i)=matmul(h_inver,v(:,i))
      enddo 
      call pot(f,frc,s,ntot,h,str)

! 温度圧力の計算
      ekin=0.d0
      do i=1,ntot
       ekin=ekin+0.5d0*mass(i)*sum(v(:,i)**2)
      enddo
      tempk=ekin*2.d0/3.d0/ntot/kb
      call p_tesor(s,ds,h,dh,str,ntot,mass,vol,p,e_kin)
      write(*,*)tempk,p(1,1)*Gp
!
      Preg=0.d0
      Preg(1,1)=1d-4/Gp; Preg(2,2)=1d-4/Gp; Preg(3,3)=1d-4/Gp 
      dh=0.d0
      w=5.d6

!=========================計算スタート====================
      do istep=1,maxstep

!　速度の更新
      GinvGdot=matmul(h_inver,dh)+transpose(matmul(h_inver,dh))
 !     G=matmul(h_inver,dh)
 !     GinvGdot_=matmul(matmul(h_inver,G),h)+G
 !     write(*,*)GinvGdot_-GinvGdot   !1d-10の差　丸め誤差？
      
      do i=1,ntot
        ds(:,i)=ds(:,i)+0.5d0*dt*matmul(h_inver,frc(:,i)/mass(i))
     &                 -0.5d0*dt*matmul(GinvGdot,ds(:,i))
      enddo
      call p_tesor(s,ds,h,dh,str,ntot,mass,vol,p,e_kin)
      dh=dh+0.5d0*dt*matmul((p-Preg),sgm)/w
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
      call vol_inverse(h,vol,h_inver,sgm)
      call pot(f,frc,s,ntot,h,str)

!　速度の更新
      GinvGdot=matmul(h_inver,dh)+transpose(matmul(h_inver,dh))
      do i=1,ntot
        ds(:,i)=ds(:,i)+0.5d0*dt*matmul(h_inver,frc(:,i)/mass(i))
     &                 -0.5d0*dt*matmul(GinvGdot,ds(:,i))
      enddo
      call p_tesor(s,ds,h,dh,str,ntot,mass,vol,p,e_kin)
      dh=dh+0.5d0*dt*matmul((p-Preg),sgm)/w
!
!　温度・圧力の計算
      ekin=(e_kin(1,1)+e_kin(2,2)+e_kin(3,3))*0.5d0
      tempk=ekin*2.d0/3.d0/ntot/kb
      gpa(istep)=(p(1,1)+p(2,2)+p(3,3))/3.d0
!

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
      write(*,*)istep,tempk,gpa(istep)*Gp,h(1,1),h(2,2),h(3,3),hamil
      
      rec1(istep)=tempk
      rec2(istep)=f
      rec3(istep)=ekin
      rec4(istep)=0.5d0*w*sum(dh**2)
      rec5(istep)=vol*Preg(1,1)
      rec6(istep)=hamil
      
      enddo

!=====================loopここまで==============================

!　グラフ
      smap=0.d0
      call sma(gpa,smap,maxstep)
      open(12,file='t.dat')
      do i=1,maxstep
        write(12,*)rec1(i),rec2(i),rec3(i),rec4(i),rec5(i),rec6(i)
     &             ,gpa(i)*gp,smap(i)*Gp
      enddo
      close(12)  
!
!　記録
      do i=1,ntot
        x(:,i)=matmul(h,s(:,i))
        v(:,i)=matmul(h,ds(:,i))+matmul(dh,s(:,i))
      enddo
 !     call kiroku(h,ntot,bohr,x,v)
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
      subroutine pot(f,frc,s,ntot,h,str)
      implicit real*8 (a-h,o-z)
      integer hx_lc,hy_lc,hz_lc,hxyz_lc,hyz_lc,pea
      real*8 frc(3,ntot),s(3,ntot),h(3,3),str(3,3),rij(3)
      parameter(sgm=3.4d0/0.5292d0)  
      parameter(eps=120d0/11605d0/27.2116d0)
      parameter(sgm12=sgm**12,sgm6=sgm**6)
      parameter(cutoff=2.5d0*sgm)
      integer,allocatable,dimension(:)::lshd,lscl

      frc=0.d0
      str=0.d0
      pea=0

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
                rij(:)=rij(:)-dnint(rij(:))
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
      subroutine p_tesor(s,ds,h,dh,str,ntot,mass,vol,p,e_kin)
      implicit none
      integer ntot,i,k,l
      real*8 s(3,ntot),ds(3,ntot),str(3,3),p(3,3),h(3,3)
      real*8 e_kin(3,3),mass(ntot),vol,vi(3),dh(3,3)
            
      e_kin=0.d0
      do i=1,ntot
        vi=matmul(h,ds(:,i))+matmul(dh,s(:,i))
       do k=1,3
       do l=1,3  
       e_kin(k,l)=e_kin(k,l)+mass(i)*vi(k)*vi(l)
       enddo
       enddo
      enddo
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