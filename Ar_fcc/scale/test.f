      program ljmd
      
      implicit real*8 (a-h,o-z)
      real*8 kb
      parameter(ev_k=11605.d0)
      parameter(hartree_ev=27.2116d0)
      parameter(kb=1.d0/ev_k/hartree_ev)
      parameter(gpa=29421.d0)
      parameter(bohr=0.5292d0)

      real*8,allocatable,dimension(:,:)::x,v,frc,s,ds
      real*8,allocatable,dimension(:)::mass
      character*8,allocatable,dimension(:)::lsp
      character*40 filename
      real*8 p(3,3),h(3,3),vh(3,3),h_inver(3,3)
      real*8 G(3,3),e_kin(3,3),str(3,3)
      integer neighbor_counter(13500),neighbor_list(300,13500)

      real*8,allocatable,dimension(:)::rec1,rec2,rec3,rec4,rec5,rec6
      
!　.datファイルの読み込み      
      open(10,file='lv50_bs1d-4.dat')
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

      call ol_inverse(h,vol,h_inver)

!　スケール  
      do i=1,ntot
        s(:,i)=matmul(h_inver,x(:,i))
        ds(:,i)=matmul(h_inver,v(:,i))
      enddo 

! 温度・圧力の計算
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
      call neighbor(s,ntot,h,neighbor_counter,neighbor_list)
      call pot(ntot,s,h,str,frc,neighbor_list,neighbor_counter)
      write(*,*)ntot,str(1,1)

      



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
!　linked cell
      subroutine neighbor
     &   (s,ntot,h,neighbor_counter,neighbor_list)
      implicit none
      integer mx,my,mz,hx_lc,hy_lc,hz_lc,m1,m,ntot
      real*8 h(3,3),s(3,ntot)
      real*8 cutoff,rij2
      integer hyz_lc,hxyz_lc,i,m1x,m1y,m1z
      integer ishiftx,ishifty,ishiftz,j
      integer kuxmax,kuymax,kuzmax,kux,kuy,kuz
      integer,allocatable,dimension(:):: lshd,lscl
      real*8 rij(3),hx_cell,hy_cell,hz_cell
      integer max_neighbor
      parameter (max_neighbor=200)
      integer neighbor_counter(ntot),neighbor_list(max_neighbor,ntot)
      parameter(cutoff=2.5d0*3.4d0/0.5292d0) !cutoff=16.06

       neighbor_counter=0

       hx_lc=ceiling(sqrt(dot_product(h(:,1),h(:,1)))/cutoff)-1 !x座標のセル数
       hy_lc=ceiling(sqrt(dot_product(h(:,2),h(:,2)))/cutoff)-1
       hz_lc=ceiling(sqrt(dot_product(h(:,3),h(:,3)))/cutoff)-1
       hyz_lc=hy_lc*hz_lc     !セルのyz平面の数
       hxyz_lc=hyz_lc*hx_lc   !総セル数 
       hx_cell=1.d0/hx_lc     !各方向のセルの長さ ボックスサイズは１
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
                if (rij2<(cutoff*1.2)**2) then
                 neighbor_counter(i)=neighbor_counter(i)+1
                 neighbor_counter(j)=neighbor_counter(j)+1
                 neighbor_list(neighbor_counter(i),i)=j
                 neighbor_list(neighbor_counter(j),j)=i
                endif
              end if
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
! call pot      
      subroutine pot(ntot,s,h,virial,dfdx,neighbor_list,neighbor_counter)
      implicit none
      integer ntot,j,k,l,a,b
      real*8 sgm,eps,sgm12,sgm6,cutoff,rij(3),rij2
      real*8 f,dfdx(3,ntot),factor,virial(3,3),s(3,ntot),h(3,3)
      integer neighbor_counter(ntot),max_neighbor
      parameter (max_neighbor=200)
      integer neighbor_list(max_neighbor,ntot)

      parameter(sgm=3.4d0/0.5292d0)
      parameter(eps=120d0/11605d0/27.2116d0)
      parameter(sgm12=sgm**12,sgm6=sgm**6)
      parameter(cutoff=2.5d0*sgm)

      dfdx=0d0
      f=0d0
      virial=0d0

      do a=1,ntot
      do b=1,neighbor_counter(a)
       j=neighbor_list(b,a)
       if (a<j)then
       rij=s(:,a)-s(:,j)
       rij=rij-dnint(rij)
       rij=matmul(h,rij)
       rij2=dot_product(rij,rij)
       if (rij2<cutoff**2)then
        f=f+4d0*eps*(sgm12/rij2**6-sgm6/rij2**3)
        factor=4d0*eps*
     &       (-12d0*sgm12/rij2**7+6d0*sgm6/rij2**4)
        dfdx(1,a)=dfdx(1,a)-factor*rij(1)
        dfdx(2,a)=dfdx(2,a)-factor*rij(2)
        dfdx(3,a)=dfdx(3,a)-factor*rij(3)
        dfdx(1,j)=dfdx(1,j)+factor*rij(1)
        dfdx(2,j)=dfdx(2,j)+factor*rij(2)
        dfdx(3,j)=dfdx(3,j)+factor*rij(3)
        do k=1,3
        do l=1,3
         virial(k,l)=virial(k,l)-factor*rij(k)*rij(l)
        enddo
        enddo
       endif
       endif
      enddo
      enddo

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
!
  
 
      




