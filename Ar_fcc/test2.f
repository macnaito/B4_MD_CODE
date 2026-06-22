      program pbc_lvBAOAB_parrinello
! PBC lcl lengevine Parrinello-Rahman
      implicit none
      integer nmax
      parameter(nmax=4000)
      integer i,j,k,m,n,maxstep,itemp,istep
      real*8 x(3,nmax),v(3,nmax),dfdx(3,nmax)
      real*8 f,dt
      real*8 h(3,3),vol
      real*8 amass,bohr
      character*2 lsp(nmax)
      character*40 filename2
      parameter(amass=40d0*1836d0)
      parameter(bohr=0.5292d0)
      real*8 W,sgm(3,3)
      real*8 p(3,3),kb,pext(3,3)
      real*8 s(3,nmax),ah(3,3),vh(3,3),h_inver(3,3),dmt(3,3)
      real*8 vs(3,nmax),as(3,nmax),mt(3,3),mt_inver(3,3)
      
! time_step      
      parameter(maxstep=20)
      real*8 gpa(maxstep)
      real*8 t(maxstep),tk,kin,rmin
    

! rang
      real*8 gamma,Treg
      real*8 sigma
      real*8 gauss
      external gauss
      dt=41d0
      gamma=1.d-3
      kb=1.d0/(27.2116*11605d0)
      pext=0d0
     

!ファイルの読みこみ      
      open(10,file='lv50k_bs0.001gpa.dat')
 !     open(10,file='init.dat')
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

      

      filename2='out000.xyz'
      k=0
      pext=0d0
    
      ah=0d0
      vh=0d0
      W=1.d10
      as=0d0


!計算start    
      call pot (f,dfdx,x,n,h) ! mt:metoric tensor
      call inverse_mass(h,h_inver,vol,mt,mt_inver,sgm)
      call press(v,dfdx,x,n,vol,p)

      do i=1,n
       s(:,i)=matmul(h_inver,x(:,i)) 
       vs(:,i)=matmul(h_inver,v(:,i))
      enddo
    
      dmt=matmul(transpose(vh),h)+matmul(transpose(h),vh)
      as=matmul(h_inver,dfdx)/amass-matmul(mt_inver,(matmul(dmt,vs)))
      ah=1.d0/W*matmul((p-pext),sgm)

      do istep=1,maxstep
        Treg=50d0
        sigma=sqrt(kb*Treg/amass*(1d0-exp(-gamma*dt)**2))

        do i=1,n
          vs(:,i)=vs(:,i)+0.5d0*dt*as(:,i)
        enddo

        vh=vh+0.5d0*dt*ah
        
   
        do i=1,n
         s(:,i)=s(:,i)+dt*vs(:,i)
        enddo

        do j=1,3
          do i=1,n
           if(s(j,i) >= 1.d0) s(j,i)=s(j,i)-1.d0
           if(s(j,i) < 0.d0)  s(j,i)=s(j,i)+1.d0
          enddo
        enddo

        h=h+dt*vh

        x=matmul(h,s)
        v=matmul(h,vs)
!        do j=1,n
!          v(:,j)=exp(-gamma*dt)*v(:,j)+sigma*gauss()
!        enddo
!        vs=matmul(h_inver,v-matmul(vh,s))

        call pot(f,dfdx,x,n,h)
        call inverse_mass(h,h_inver,vol,mt,mt_inver,sgm)
        call press(v,dfdx,x,n,vol,p)


      dmt=matmul(transpose(vh),h)+matmul(transpose(h),vh)
      as=matmul(h_inver,dfdx)/amass-matmul(mt_inver,(matmul(dmt,vs)))
      ah=1.d0/W*matmul((p-pext),sgm)

      do i=1,n
       vs(:,i)=vs(:,i)+0.5d0*dt*as(:,i)
      enddo
      vh=vh+0.5d0*dt*ah      
     
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
            write(11,'(a2,i5,4e15.7)')lsp(m),m,
     &       x(1,m)*bohr,x(2,m)*bohr,x(3,m)*bohr,0d0
          enddo
          close(11)
        endif
!ここまで      

!温度計算    
        v=matmul(h,vs)+matmul(vh,s)    
        kin=0d0
        do j=1,n
          kin=kin+0.5d0*amass*
     &     (v(1,j)**2+v(2,j)**2+v(3,j)**2)
        enddo
        tk=kin*2d0/(3d0*dble(n))*
     &     27.2116*11605d0
!温度計算 
 
        write(*,*) istep,h(1,1),tk,ah(1,1),rmin
        gpa(istep)=p(1,1)*29421d0
        t(istep)=tk
        t(istep)=tk
      enddo

!グラフ      
      open (12,file='t.dat')
      do i=1,maxstep
        write(12,*) gpa(i),t(i)
      enddo
      close(12)
!ここまで

      end program

! call pot      
      subroutine pot(f,dfdx,x,n,h)
      implicit none

      integer n,i,j
      real*8 f,x(3,n),dfdx(3,n)
      real*8 xij,yij,zij,r2,factor
      real*8 sgm,eps,sgm12,sgm6
      real*8 h(3,3),cutoff
      parameter(sgm=3.4d0/0.5292d0)
      parameter(eps=120d0/11605d0/27.2116d0)
      parameter(sgm12=sgm**12,sgm6=sgm**6)
      parameter(cutoff=2.5d0*sgm)
      f=0d0
      dfdx=0d0
c-----potential
      do i=1,n-1
        do j=i+1,n
          xij=x(1,i)-x(1,j)
          yij=x(2,i)-x(2,j)
          zij=x(3,i)-x(3,j)
          xij=xij-h(1,1)*dnint(xij/h(1,1))
          yij=yij-h(2,2)*dnint(yij/h(2,2))
          zij=zij-h(3,3)*dnint(zij/h(3,3))
          r2=xij**2+yij**2+zij**2
          if (r2 < cutoff**2) then
            f=f+4d0*eps*(sgm12/r2**6-sgm6/r2**3)
          endif
        enddo
      enddo

c-----force
      do i=1,n
        do j=1,n
          if(j.ne.i)then
            xij=x(1,i)-x(1,j)
            yij=x(2,i)-x(2,j)
            zij=x(3,i)-x(3,j)
            xij=xij-h(1,1)*dnint(xij/h(1,1))
            yij=yij-h(2,2)*dnint(yij/h(2,2))
            zij=zij-h(3,3)*dnint(zij/h(3,3))
            r2=xij**2+yij**2+zij**2
            if (r2 < cutoff**2) then
              factor=4d0*eps*
     &        (-12d0*sgm12/r2**7+6d0*sgm6/r2**4)
              dfdx(1,i)=dfdx(1,i)-factor*xij
              dfdx(2,i)=dfdx(2,i)-factor*yij
              dfdx(3,i)=dfdx(3,i)-factor*zij
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
      subroutine press(v,dfdx,x,n,vol,p)
      implicit none
      real*8 v(3,n),p(3,3),dfdx(3,n),x(3,n)
      real*8 vol
      real*8 amass,a(3,3)
      integer n
      parameter(amass=40d0*1836d0)

      a=amass*matmul(v,transpose(v))
      p=1.d0/vol*(a+matmul(dfdx,transpose(x)))

      return
      end subroutine
