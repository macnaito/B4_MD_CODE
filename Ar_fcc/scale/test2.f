!module 
      module ewald_module_s
 !subroutine-program間で共有する変数
      implicit none
      real*8 bohr
      real*8 h_inver(3,3),sigm(3,3)
      parameter(pi=3.141592653589793d0)
      integer lc(5),kumax(3)
      parameter(nbin=10000)
      parameter(rc_real=10.d0)
      parameter(bohr=0.52918d0)
      real*8,allocatable,dimension(:,:)::vs,ak
      real*8,allocatable,dimension(:)::uk
      integer,allocatable,dimension(:)::lshd,lscl
      
      contains
 !

 !subroutine first
      subroutine first(h)
      implicit none
      real*8 accurary,dr2,r,erfv,h(3,3)
      real*8 b(3,3),gnorm(3),gmin
      integer i
      parameter(accurary=1d-5)

      real*8 ak2
      integer ltmp,lmax(3),knum,lx,ly,lz


      h_(1)=sqrt(dot_product(h(:,1),h(:,1)))
      h_(2)=sqrt(dot_product(h(:,2),h(:,2)))
      h_(3)=sqrt(dot_product(h(:,3),h(:,3)))
      b=2.d0*pi*transpose(h_inver)
      do i=1,3
       gnorm(i)=sqrt(dot_product(b(:,i),b(:,i)))
      enddo
      gmin=min(gnorm(1),min(gnorm(2),gnorm(3)))

      gamma=rc_real/(2.d0*(sqrt(-log(accurary))))
      dr2=rc_real**2/nbin
      dr2i=1.d0/dr2

  !実空間のポテンシャルと力のテーブル
      if(allocated(ak)) deallocate(ak,uk,vs)
      allocate(vs(nbin,0:1))
      do i=1,nbin
       r=sqrt(dble(i)*dr2)
       erfv=erfc(r*0.5d0/gamma)
       vs(i,0)=erfv/r         
       vs(i,1)=erfv/r**3+
     &    1.d0/sqrt(pi)/gamma/r**2*exp(-(r*0.5d0/gamma)**2)
      enddo

      do i=1,nbin   !curoff補正
       r=sqrt(dble(i)*dr2) 
       vs(i,0)=vs(i,0)-vs(nbin,0)+(r-rc_real)*rc_real*vs(nbin,1)
       vs(i,1)=vs(i,1)-rc_real/r*vs(nbin,1)
      enddo
 
  !kベクトル　uk
      do i=1,3
       ltmp=1
       do while(exp(-(gamma*gmin*ltmp)**2)/ltmp**2>accurary)
        ltmp=ltmp+1
       enddo
       lmax(i)=ltmp
      enddo

      knum=0
      do lx=-lmax(1),lmax(1)
      do ly=-lmax(2),lmax(2)
      do lz=0,lmax(3)
        ak2=dot_product(lx*b(:,1)+ly*b(:,2)+lz*b(:,3),
     &                  lx*b(:,1)+ly*b(:,2)+lz*b(:,3))
        if(ak2>0.d0) then
          if(exp(-gamma**2*ak2)/ak2>accurary) then
            knum=knum+1
          endif
        endif
      enddo
      enddo
      enddo
      kmax=knum

      allocate(ak(3,kmax),uk(kmax))
      knum=0
      do lx=-lmax(1),lmax(1)
      do ly=-lmax(2),lmax(2)
      do lz=0,lmax(3)
        ak2=dot_product(lx*b(:,1)+ly*b(:,2)+lz*b(:,3),
     &                  lx*b(:,1)+ly*b(:,2)+lz*b(:,3))
        if(ak2>0.d0) then
          if(exp(-gamma**2*ak2)/ak2>accurary) then
          knum=knum+1
          ak(:,knum)=lx*b(:,1)+ly*b(:,2)+lz*b(:,3)
          uk(knum)=4.d0*pi/vol/ak2*exp(-gamma**2*ak2)
          if(lz>0) uk(knum)=2.d0*uk(knum)
          endif
        endif
      enddo
      enddo
      enddo

      return
      end subroutine
 !

 !subroutine linked-cell-list
      subroutine lcl(ntot,s)      
      implicit none
      integer i,ntot,mm
      integer m(3)
      real*8 rc(3)
      real*8 s(3,ntot)

      lc(1:3)=max(int(h_(1:3)/rc_real),1)
      lc(4)=lc(2)*lc(3)
      lc(5)=lc(1)*lc(4)  !総セル数
      rc(1:3)=1.d0/lc(1:3)   !セルの長さ
      allocate(lshd(lc(5)),lscl(ntot))
      lshd=0

      do i=1,ntot
       m(1:3)=int(s(1:3,i)/rc(1:3))
       m(1:3)=min(max(m(1:3),0),lc(1:3)-1)
       mm=m(1)*lc(4)+m(2)*lc(3)+m(3)+1
       lscl(i)=lshd(mm)
       lshd(mm)=i
      enddo

      return
      end
 !

 !subroutine real-space froces energy stress
      subroutine real(ntot,s,h,qi,lsp,kb,frc_r,epot_r,str_r,rmin,
     &              epot_lj,str_lj,frc_lj)
      implicit none
      integer ntot,i,j,ir,ir1,mx,my,mz,m,m1,np,u,v
      integer kux,kuy,kuz,m1x,m1y,m1z,ishiftx,ishifty,ishiftz
      real*8 rij(3),rij2,s(3,ntot),fr,rmin,kb,h(3,3)
      real*8 qi(ntot),vs0,epot_r,vs1,es1,frc_r(3,ntot),st1,str_r(3,3)
      real*8 sgm,sgm_na,sgm_cl,eps,eps_na,eps_cl,sgm6,sgm12
      real*8 epot_lj,factor,frc_lj(3,ntot),str_lj(3,3)
      character*2 lsp(ntot)
      !パラメータ joung-cheathamモデル
      parameter(sgm_na=2.4392d0)   !単位Å
      parameter(sgm_cl=4.4172d0)
      parameter(eps_na=49.98d0)     !K
      parameter(eps_cl=50.32d0)

      epot_r=0.d0
      frc_r=0.d0
      str_r=0.d0
      epot_lj=0.d0
      frc_lj=0.d0
      str_lj=0.d0
      rmin=100.d0
      np=0

      kumax(:)=int(rc_real/h_(:)+1)
      do mz=0,lc(3)-1
      do my=0,lc(2)-1
      do mx=0,lc(1)-1
       m=mx*lc(4)+my*lc(3)+mz+1
       if(lshd(m)==0) cycle
       do kuz=-kumax(3),kumax(3)
       do kuy=-kumax(2),kumax(2)
       do kux=-kumax(1),kumax(1)
        m1x=mx+kux
        m1y=my+kuy
        m1z=mz+kuz
        call get_ishift_ewald(lc(1),m1x,ishiftx)
        call get_ishift_ewald(lc(2),m1y,ishifty)
        call get_ishift_ewald(lc(3),m1z,ishiftz)
        m1=m1x*lc(4)+m1y*lc(3)+m1z+1
        if(lshd(m1)==0) cycle
         i=lshd(m)
         do while(i>0)
          j=lshd(m1)
          do while(j>0)
          if(i<j) then
            rij(:)=s(:,i)-s(:,j)
            rij(:)=rij(:)-dnint(rij(:))
            rij(:)=matmul(h,rij(:))
            rij2=sum(rij(:)**2)
            rmin=min(rmin,rij2)
            if(rij2<rc_real**2)then
 !ここからはewald法
             np=np+1
             ir=int(rij2*dr2i)
             if(ir==0) then 
             fr=0.d0
             else
             fr=(rij2*dr2i)-ir
             endif
             ir=max(1,ir)
             ir1=min(ir+1,nbin) 
             vs0=(1.d0-fr)*vs(ir,0)+fr*vs(ir1,0)
             epot_r=epot_r+qi(i)*qi(j)*vs0
             vs1=(1.d0-fr)*vs(ir,1)+fr*vs(ir1,1)
             es1=qi(i)*qi(j)*vs1
             frc_r(:,i)=frc_r(:,i)+rij(:)*es1
             frc_r(:,j)=frc_r(:,j)-rij(:)*es1
             st1=es1/vol
             str_r(:,1)=str_r(:,1)+rij(:)*rij(1)*st1
             str_r(:,2)=str_r(:,2)+rij(:)*rij(2)*st1
             str_r(:,3)=str_r(:,3)+rij(:)*rij(3)*st1
 !ここからはLJ
             if(lsp(i)=='Na'.and.lsp(j)=='Na') then
              sgm=sgm_na/bohr
              eps=eps_na*kb
             elseif(lsp(i)=='Cl'.and.lsp(j)=='Cl') then
              sgm=sgm_cl/bohr
              eps=eps_cl*kb
             else
               sgm=(sgm_na+sgm_cl)/bohr/2.d0
               eps=sqrt(eps_na*eps_cl)*kb
             endif
             sgm12=sgm**12
             sgm6=sgm**6
             epot_lj=epot_lj+4.d0*eps*(sgm12/(rij2**6)-sgm6/(rij2**3))
             factor=4.d0*eps*((-12.d0*sgm12)/(rij2**7)
     &                    +(6.d0*sgm6)/(rij2**4))
             frc_lj(:,i)=frc_lj(:,i)-factor*rij(:)
             frc_lj(:,j)=frc_lj(:,j)+factor*rij(:)
             do u=1,3
             do v=1,3
              str_lj(u,v)=str_lj(u,v)-factor*rij(U)*rij(v)/vol
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
   !   write(*,*)np
      deallocate(lshd,lscl)

      return
      end
 !

 !subroutine wave-number-space forces energy stress
      subroutine wave_number(ntot,s,h,qi,frc_w,epot_w,str_w)
      implicit none
      integer ntot,k,i,m,n
      real*8 qk(2,ntot),ftemp(3,ntot),s(3,ntot),qi(ntot),h(3,3)
      real*8 csum,ssum,fkx,fky,sumvk,scv,ak2_,arg
      real*8 epot_w,frc_w(3,ntot),str_w(3,3)
      
      epot_w=0.d0
      frc_w=0.d0
      str_w=0.d0
      

      ftemp=0.d0
      do k=1,kmax
       csum=0.d0
       ssum=0.d0
       do i=1,ntot
        arg=dot_product(ak(:,k),matmul(h,s(:,i)))
        fkx=qi(i)*cos(arg)
        fky=qi(i)*sin(arg)
        qk(1,i)=fkx
        qk(2,i)=fky
        csum=csum+fkx
        ssum=ssum+fky
       enddo
       sumvk=uk(k)*(csum**2+ssum**2)
       epot_w=epot_w+0.5d0*sumvk
       do i=1,ntot
        scv=uk(k)*(qk(2,i)*csum-qk(1,i)*ssum)
        ftemp(:,i)=ftemp(:,i)+scv*ak(:,k)
       enddo
       sumvk=sumvk/vol
       ak2_=1.d0/(ak(1,k)**2+ak(2,k)**2+ak(3,k)**2)
       do m=1,3
       do n=1,3
        str_w(m,n)=str_w(m,n)+
     &   sumvk*(0.5d0-(ak2_+gamma**2)*ak(m,k)*ak(n,k))
       enddo
       enddo
      enddo 
      do i=1,ntot
       frc_w(:,i)=frc_w(:,i)+ftemp(:,i)
      enddo
      
      return
      end subroutine
 !

 !subroutine self-terms forces energy stress
      subroutine self_terms(ntot,qi,epot_s,frc_s,str_s)
      implicit none

      integer ntot,i,kux,kuy,kuz,ir,ir1,m,n
      real*8 epot_s,frc_s(3,ntot),str_s(3,3),qi(ntot)
      real*8 q2sum,rij(3),rij2,fr,vs0,vs1,st1

      epot_s=0.d0
      frc_s=0.d0
      str_s=0.d0

      do i=1,ntot
       epot_s=epot_s-0.5d0*qi(i)**2/sqrt(pi)/gamma
      enddo
      q2sum=0.d0
      do i=1,ntot
       q2sum=q2sum+qi(i)**2
      enddo
      do kux=-kumax(1),kumax(1)
      do kuy=-kumax(2),kumax(2)
      do kuz=-kumax(3),kumax(3)
       if(kux==0.and.kuy==0.and.kuz==0) cycle
       rij(:)=[kux*h_(1),kuy*h_(2),kuz*h_(3)]
       rij2=sum(rij(:)**2)
       if(rij2<rc_real**2) then
        ir=int(rij2*dr2i)
        if(ir==0) then 
         fr=0.d0
        else
         fr=(rij2*dr2i)-ir
        endif
        ir=max(1,ir)
        ir1=min(ir+1,nbin)
        vs0=(1.d0-fr)*vs(ir,0)+fr*vs(ir1,0)
        epot_s=epot_s+0.5d0*q2sum*vs0
        vs1=(1.d0-fr)*vs(ir,1)+fr*vs(ir1,1)
        st1=vs1*0.5d0*q2sum/vol
        do m=1,3
        do n=1,3
         str_s(m,n)=str_s(m,n)+rij(m)*rij(n)*st1
        enddo
        enddo
       endif
      enddo
      enddo
      enddo 

      return
      end
 !

 !subroutine shift
      subroutine get_ishift_ewald(lc,i,ishift)
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

 ! maxwell分布
      subroutine gauss(ga)
      implicit none
      real*8 ga,u1,u2

      call random_number(u1)
      call random_number(u2)
      if(u1<1d-12) u1=1d-12
      ga=sqrt(-2.d0*log(u1))*cos(2.d0*3.14159265d0*u2)
      return
      end
 !

 ! 逆行列と体積
      subroutine inverse_and_mass(h)
      implicit real*8 (a-h,o-z)
      real*8 h(3,3)
      integer i

      vol=h(1,1)*(h(2,2)*h(3,3)-h(2,3)*h(3,2))
     &   -h(1,2)*(h(2,1)*h(3,3)-h(2,3)*h(3,1))
     &   +h(1,3)*(h(2,1)*h(3,2)-h(2,2)*h(3,1))
      do i=1,3
       h_inver(1,i)=(h(2,mod(i,3)+1)*h(3,mod(i+1,3)+1)
     &              -h(2,mod(i+1,3)+1)*h(3,mod(i,3)+1))/vol
       h_inver(2,i)=(h(1,mod(i+1,3)+1)*h(3,mod(i,3)+1)
     &              -h(1,mod(i,3)+1)*h(3,mod(i+1,3)+1))/vol
       h_inver(3,i)=(h(1,mod(i,3)+1)*h(2,mod(i+1,3)+1)
     &              -h(1,mod(i+1,3)+1)*h(2,mod(i,3)+1))/vol
      enddo 

      sigm(1,1)=(h(2,2)*h(3,3)-h(2,3)*h(3,2))
      sigm(1,2)=(h(1,3)*h(3,2)-h(1,2)*h(3,3))
      sigm(1,3)=(h(1,2)*h(2,3)-h(1,3)*h(2,2))
      sigm(2,1)=(h(2,3)*h(3,1)-h(2,1)*h(3,3))
      sigm(2,2)=(h(1,1)*h(3,3)-h(1,3)*h(3,1))
      sigm(2,3)=(h(1,3)*h(2,1)-h(1,1)*h(2,3))
      sigm(3,1)=(h(2,1)*h(3,2)-h(2,2)*h(3,1))
      sigm(3,2)=(h(1,2)*h(3,1)-h(1,1)*h(3,2))
      sigm(3,3)=(h(1,1)*h(2,2)-h(1,2)*h(2,1))
     
      return
      end
 !

      end module ewald_module_s
!

      
!?!?!--------------ここからmain program--------------!?!?!
      program ewald
      use ewald_module_s
!  NaClのMD lj-ewald法 linked-cell-list 一般化座標　エネルギー保存
! エネルギーは6d-3変化
      implicit real*8 (a-h,o-z)
      parameter(eV_K=11605.d0) !K
      parameter(Hartree_eV=27.2116d0) !eV
      parameter(Gpa=29421.d0)

      real*8 str(3,3),str_r(3,3),str_w(3,3),str_s(3,3),h(3,3)
      real*8 str_lj(3,3)
      real*8,allocatable,dimension(:,:)::x,v,frc,frc_r,frc_w,
     &                                   frc_s,frc_lj,s,ds
      real*8,allocatable,dimension(:)::qi,mass
      character*2,allocatable,dimension(:)::lsp
      character*40 filename
      real*8 kinetic_energy,kb,G(3,3)
      real*8 ekin(3,3),p(3,3)
    
      
      integer maxstep
      parameter(maxstep=2000)
      real*8 rec1(maxstep),rec2(maxstep),rec3(maxstep)
      real*8 rec4(maxstep),rec5(maxstep),rec6(maxstep)

      dt=41.34d0*0.1d0  !単位 au　1step:0.005fs
      k=0
      
      call cpu_time(t1)

!ファイルの読み込み
      open(10,file='init.dat')
      read(10,*)h(1,1),h(2,1),h(3,1) !単位Å
      read(10,*)h(1,2),h(2,2),h(3,2)
      read(10,*)h(1,3),h(2,3),h(3,3)
      read(10,*)ntot
       allocate(x(3,ntot),v(3,ntot))
       allocate(qi(ntot),lsp(ntot),mass(ntot))
      do i=1,ntot     !単位はÅ
        read(10,*)lsp(i),qi(i),x(1,i),x(2,i),x(3,i)
      enddo
      read(10,*)h(1,1),h(2,1),h(3,1) !単位Å
      read(10,*)h(1,2),h(2,2),h(3,2)
      read(10,*)h(1,3),h(2,3),h(3,3)
      close(10)

      allocate(frc(3,ntot),s(3,ntot),ds(3,ntot))
      h=h/bohr
      x(:,:)=x(:,:)/bohr   !auに変換
      do i=1,ntot   !質量の設定　単位 au
        if(lsp(i)=='Na') mass(i)=23d0*1836d0
        if(lsp(i)=='Cl') mass(i)=35.5*1836d0
        if(lsp(i)=='Ar') mass(i)=40d0*1836d0
      enddo
!

!初速度
      treg=100.d0
      kb=1.d0/Hartree_eV/eV_K
      do i=1,ntot
        do j=1,3
         call gauss(ga)
         v(j,i)=sqrt(kb*treg/mass(i))*ga
        enddo
      enddo

      kinetic_energy=0.d0
      do i=1,ntot
       kinetic_energy=kinetic_energy+0.5d0*mass(i)*sum(v(:,i)**2)
      enddo
      tempK=kinetic_energy*2.d0/(3.d0*ntot)*Hartree_eV*eV_K
      write(*,*)tempK,kinetic_energy
!

!一般化座標
      call inverse_and_mass(h)
      do i=1,ntot
       s(:,i)=matmul(h_inver,x(:,i))
       ds(:,i)=matmul(h_inver,v(:,i))
      enddo


      call first(h)
      call lcl(ntot,s)
      call real(ntot,s,h,qi,lsp,kb,frc_r,epot_r,str_r,rmin,
     &              epot_lj,str_lj,frc_lj)
      call wave_number(ntot,s,h,qi,frc_w,epot_w,str_w)
      call self_terms(ntot,qi,epot_s,frc_s,str_s)
      frc=frc_r+frc_w+frc_s+frc_lj
      epot=epot_r+epot_w+epot_s+epot_lj
      str=str_r+str_w+str_s+str_lj

 !圧力計算
      ekin=0.d0
      do i=1,ntot
       do l=1,3
       do m=1,3
        ekin(l,m)=ekin(l,m)+mass(i)*v(l,i)*v(m,i)
       enddo
       enddo 
      enddo
      do l=1,3
      do m=1,3
        p(l,m)=ekin(l,m)/vol+str(l,m)
      enddo
      enddo
      write(*,*)(p(1,1)+p(2,2)+p(3,3))*Gpa,'gpa'
 !
 !    
          


!------------------------計算start-----------------------
      do istep=1,maxstep

        G=matmul(transpose(h),h)
        kinetic_energy=0.d0
        do i=1,ntot
         kinetic_energy=kinetic_energy+0.5d0*mass(i)
     &    *dot_product(ds(:,i),matmul(G,ds(:,i)))
        enddo
        xi=xi+0.5d0*dt*(2.d0*kinetic_energy-gkbt)/Q 

 ! 1/2速度更新   
        v=matmul(h,ds)
        do i=1,ntot
         v(:,i)=v(:,i)+0.5d0*dt*frc(:,i)/mass(i)
         ds(:,i)=matmul(h_inver,v(:,i))
        enddo
 !

 ! 位置更新＋pbc        
        do i=1,ntot
         s(:,i)=s(:,i)+dt*ds(:,i)
         do j=1,3
          if(s(j,i)>=1.d0) s(j,i)=s(j,i)-1.d0
          if(s(j,i)< 0.d0) s(j,i)=s(j,i)+1.d0
         enddo
        enddo 
 ! 
        call inverse_and_mass(h)
        call lcl(ntot,s)
        call real(ntot,s,h,qi,lsp,kb,frc_r,epot_r,str_r,rmin,
     &              epot_lj,str_lj,frc_lj)
        call wave_number(ntot,s,h,qi,frc_w,epot_w,str_w)
        call self_terms(ntot,qi,epot_s,frc_s,str_s)
        frc=frc_r+frc_w+frc_s+frc_lj
        epot=epot_r+epot_w+epot_s+epot_lj
        str=str_r+str_w+str_s+str_lj

 !圧力計算
        ekin=0.d0
        do i=1,ntot
         do l=1,3
         do m=1,3
          ekin(l,m)=ekin(l,m)+mass(i)*v(l,i)*v(m,i)
         enddo
         enddo 
        enddo
        do l=1,3
        do m=1,3
           p(l,m)=ekin(l,m)/vol+str(l,m)
        enddo
        enddo
 !       

 ! 1/2速度更新   
        v=matmul(h,ds) 
        do i=1,ntot
         v(:,i)=v(:,i)+0.5d0*dt*frc(:,i)/mass(i)
         ds(:,i)=matmul(h_inver,v(:,i))
        enddo
 !

        kinetic_energy=0.d0
        do i=1,ntot
         kinetic_energy=kinetic_energy+0.5d0*mass(i)*sum(v(:,i)**2)
        enddo
        tempK=kinetic_energy*2.d0/(3.d0*ntot)*Hartree_eV*eV_K

        rec1(istep)=tempK
        rec2(istep)=p(1,1)*Gpa
        rec3(istep)=p(2,2)*Gpa
        rec4(istep)=p(3,3)*Gpa
        rec5(istep)=rmin
        rec6(istep)=epot+kinetic_energy

        write(*,*)istep,tempK,p(1,1)*Gpa,sqrt(rmin),epot+kinetic_energy

 ! xyzファイルへの書き出し
        filename='out000.xyz'
        if(mod(istep,int(maxstep/100))==0) then
         x=matmul(h,s) 
         k=k+1
         write(filename(4:6),'(i3.3)')k
         open(20,file=filename)
         write(20,*)ntot
         write(20,'(a,3(3e15.7),a,a)')
     &  'Lattice="',h(1,1)*bohr,h(1,2)*bohr,h(1,3)*bohr,
     &              h(2,1)*bohr,h(2,2)*bohr,h(2,3)*bohr,
     &              h(3,1)*bohr,h(3,2)*bohr,h(3,3)*bohr,'" ',
     &  'Properties=species:S:1:id:I:1:pos:R:3:qi:R:1'  
         do i=1,ntot
          write(20,'(a2,i5,3e15.7,f8.3)') 
     &          lsp(i),i,x(1,i)*bohr,x(2,i)*bohr,x(3,i)*bohr,qi(i)
         enddo
        close(20)
       endif
 !

      enddo

 !グラフ
      open(30,file='e.dat')
      do i=1,maxstep
       write(30,*)rec1(i),rec2(i),rec3(i),rec4(i),rec5(i),rec6(i)
      enddo
      close(30)
 !

      call cpu_time(t2)
      write(*,*)t2-t1

      end program ewald