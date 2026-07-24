!　スケールした座標でのMD 温度・圧力制御なし　エネルギー保存を確認
!　速度更新はvで行う。→s、dsでの更新に変更する
! リンクリストのところをsで行うように変更する
!　温度制御や圧力制御を入れる。
!
! 一種類の原子Ar だけで行う場合はmassのところを変えてもいい　今は中途半端になってる
      program ljmd
      
      implicit real*8 (a-h,o-z)
      parameter(ev_k=11605.d0)     !温度の単位変換
      parameter(hartree_ev=27.2116d0)
      real*8 kb
      parameter(kb=1.d0/ev_k/hartree_ev)
      parameter(Gp=29421.d0)       !圧力の単位変換
      parameter(bohr=0.5292d0)

      !s:scaleした座標　ds:scaleした速度dは微分の意味
      !x,v,s,dsは3*ntotの行列で表す
      real*8,allocatable,dimension(:,:)::x,v,frc,s,ds
      real*8,allocatable,dimension(:)::mass
      character*8,allocatable,dimension(:)::lsp
      character*40 filename
      real*8 p(3,3),h(3,3),h_inver(3,3)   !h_inver:逆行列
      real*8 G(3,3),e_kin(3,3),str(3,3)   !kinetic_energy運動エネルギー
                                          !e_kin:3*3の行列　ekin:スカラー

      !最後にグラフにするためのもの
      real*8,allocatable,dimension(:)::rec1,rec2,rec3,rec4,rec5,rec6
      real*8,allocatable,dimension(:)::gpa,smap

      maxstep=3000
      allocate(gpa(maxstep),smap(maxstep),rec1(maxstep),rec2(maxstep)
     &        ,rec3(maxstep),rec4(maxstep),rec5(maxstep),rec6(maxstep))

      k=0

!　.datファイルの読み込み      
      open(10,file='init.dat')   
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
!
      dt=41*5.d0

!　スケール 
      call vol_inverse(h,vol,h_inver) 
      do i=1,ntot
        s(:,i)=matmul(h_inver,x(:,i))
        ds(:,i)=matmul(h_inver,v(:,i))
      enddo 

! 温度の計算
      !いろんな方法で、温度を計算する。計算が一致するかの確認と、計算スピードを見る。
      !いらないやつは消して
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
      call pot(f,frc,s,ntot,h,str)

!=========================計算スタート====================
      do istep=1,maxstep

!　速度の更新
      v=matmul(h,ds)
      do i=1,ntot
        v(:,i)=v(:,i)+0.5d0*dt*frc(:,i)/mass(i)
        ds(:,i)=matmul(h_inver,v(:,i))
      enddo
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
      call pot(f,frc,s,ntot,h,str)

!　速度の更新
      v=matmul(h,ds)
      do i=1,ntot
        v(:,i)=v(:,i)+0.5d0*dt*frc(:,i)/mass(i)
        ds(:,i)=matmul(h_inver,v(:,i))
      enddo
! 
!　温度の計算
      ekin=0.d0
      do i=1,ntot
       ekin=ekin+0.5d0*mass(i)
     &  *(v(1,i)**2+v(2,i)**2+v(3,i)**2)
      enddo
      tempk=ekin*2.d0/3.d0/ntot/kb
!
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

      write(*,*)istep,tempk,gpa(istep),f+ekin
      
      rec1(istep)=tempk
      rec2(istep)=f
      rec3(istep)=ekin
      rec4(istep)=f+ekin
      rec5(istep)=0.d0
      rec6(istep)=0.d0
      
      enddo

!=====================loopここまで==============================

!　グラフ
      open(12,file='t.dat')
      do i=1,maxstep
        write(12,*)rec1(i),rec2(i),rec3(i),rec4(i),rec5(i),rec6(i)
     &             ,gpa(i),smap(i)
      enddo
      close(12)
      !gnuplotでusing 1~6で自分で設定したrecxの変化、using 7で圧力をみれる。  
!
!　記録
      !最後の状態（xとvをåで保存する。）
 !     call kiroku(h,ntot,bohr,x,v)
!

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

      !！！！ここをスケール座標に変える！！！！！
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

      !ここも変える！！！！！
      kuxmax=1
      kuymax=1
      kuzmax=1
      f=0d0
      
      do mz=0,hz_lc-1
      do my=0,hy_lc-1
      do mx=0,hx_lc-1
        m=mx*hyz_lc+my*hz_lc+mz+1 
        if(lshd(m)==0) cycle 
        do kuz=-kuzmax,kuzmax
        do kuy=-kuymax,kuymax
        do kux=-kuxmax,kuxmax
          m1x=mx+kux
          call get_ishift(hx_lc,m1x,ishiftx)
          m1y=my+kuy
          call get_ishift(hy_lc,m1y,ishifty)
          m1z=mz+kuz
          call get_ishift(hz_lc,m1z,ishiftz)
          m1=m1x*hyz_lc+m1y*hz_lc+m1z+1 
          if (lshd(m1)==0) cycle 
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
!　瞬間の圧力と平均の圧力
      !simple moving average（単純移動平均）
      !圧力は時間変化するため、過去nstepの平均をとる。あんまり意味ないかも。
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


