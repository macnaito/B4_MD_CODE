      program pbc_parrinello
! PBC  Parrinello-Rahman
      implicit none
      integer nmax,pair
      parameter(nmax=4000)
      integer i,j,k,m,n,maxstep,itemp,istep
      real*8 x(3,nmax),v(3,nmax),dfdx(3,nmax)
      real*8 f,ekin(3,3),dt
      real*8 h(3,3),vol
      real*8 amass,bohr
      character*2 lsp(nmax)
      character*40 filename2
      parameter(amass=40d0*1836d0)
      parameter(bohr=0.5292d0)
      real*8 W,G(3,3)
      real*8 virial(3,3),p(3,3),kb,pext(3,3)
      real*8 s(3,nmax),ah(3,3),vh(3,3),h_inver(3,3)
      real*8 vs(3,nmax),as(3,nmax),rmin
! time_step      
      parameter(maxstep=1000)
    

! rang
      real*8 gamma,Treg
      real*8 sigma
      real*8 gauss
      external gauss
      dt=41d0*5
      gamma=1.d-4
      kb=1.d0/(27.2116*11605d0)
      h=0d0

!ファイルの読みこみ      
 !     open(10,file='lv50k_bs0.001gpa.dat')
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
      read(10,*)h(1,1),h(2,1),h(3,1)
      read(10,*)h(1,2),h(2,2),h(3,2)
      read(10,*)h(1,3),h(2,3),h(3,3)
      h=h/bohr
      close(10)
!ファイルの読み込み

      filename2='out000.xyz'
      k=0
      pext=0d0
      pext(1,1)=3.399d-5
      pext(2,2)=3.399d-5
      pext(3,3)=3.399d-5
      ah=0d0
      vh=0d0
      W=1.d6
      as=0d0


!計算start    
      call pot (f,dfdx,x,n,h,virial,pair)
      call inverse_mass(h,h_inver,vol)
      call press(v,virial,n,vol,p,ekin)
 
      do i=1,n
       vs(:,i)=matmul(h_inver,v(:,i))
       s(:,i)=matmul(h_inver,x(:,i))
      enddo
      G=matmul(h_inver,vh)
      as=-matmul(h_inver,dfdx)/amass-2.d0*matmul(G,vs)
      ah=vol/W*matmul((p-pext),transpose(h_inver))


      do istep=1,maxstep
        Treg=50d0
        sigma=sqrt(kb*Treg/amass*(1d0-exp(-gamma*dt)**2))

        do i=1,n
         vs(1,i)=vs(1,i)+0.5d0*dt*as(1,i)
         vs(2,i)=vs(2,i)+0.5d0*dt*as(2,i)
         vs(3,i)=vs(3,i)+0.5d0*dt*as(3,i)
        enddo
        vh=vh+0.5d0*dt*ah
   
        do i=1,n
         s(1,i)=s(1,i)+dt*vs(1,i)
         s(2,i)=s(2,i)+dt*vs(2,i)
         s(3,i)=s(3,i)+dt*vs(3,i)
        enddo

        do j=1,3
          do i=1,n
           if(s(j,i) >= 1.d0) s(j,i)=s(j,i)-1.d0
           if(s(j,i) < 0.d0)  s(j,i)=s(j,i)+1.d0
          enddo
        enddo

        h=h+dt*vh

        x=matmul(h,s)
        v=matmul(h,vs)+matmul(vh,s)
        call pot(f,dfdx,x,n,h,virial,pair)
        call inverse_mass(h,h_inver,vol)
        call press(v,virial,n,vol,p,ekin)

        G=matmul(h_inver,vh)
        as=-matmul(h_inver,dfdx)/amass-2d0*matmul(G,vs)
        ah=vol/W*matmul((p-pext),transpose(h_inver))

        do i=1,n
         vs(1,i)=vs(1,i)+0.5d0*dt*as(1,i)
         vs(2,i)=vs(2,i)+0.5d0*dt*as(2,i)
         vs(3,i)=vs(3,i)+0.5d0*dt*as(3,i)
        enddo
        vh=vh+0.5d0*dt*ah
     
!ファイルの書き出し
        if(mod(istep,10).eq.0)then
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

        write(*,*)h(1,1),f,virial(1,1),maxval(abs(dfdx))
      enddo
      end program


!lj
      subroutine pot(f,dfdx,x,n,h,virial,pair)
      implicit none

      integer n,i,j,m,l,pair
      real*8 f,x(3,n),dfdx(3,n),virial(3,3)
      real*8 r2,factor,rij(3),h(3,3)
      real*8 sgm,eps,sgm12,sgm6,bohr,cutoff,rmin
 !-----atomic unit is used
      parameter(sgm=3.4d0/0.5292d0)
      parameter(eps=120d0/11605d0/27.2116d0)
      parameter(sgm12=sgm**12,sgm6=sgm**6)
      parameter(cutoff=2.5d0*sgm)
      parameter(bohr=0.5292d0)

      f=0d0
      virial=0d0
      rmin=10
      pair=0d0

c-----potential
      do i=1,n-1
        do j=i+1,n
          rij(1)=x(1,i)-x(1,j)
          rij(2)=x(2,i)-x(2,j)
          rij(3)=x(3,i)-x(3,j)
          rij(1)=rij(1)-h(1,1)*dnint(rij(1)/h(1,1))
          rij(2)=rij(2)-h(2,2)*dnint(rij(2)/h(2,2))
          rij(3)=rij(3)-h(3,3)*dnint(rij(3)/h(3,3))
          r2=rij(1)**2+rij(2)**2+rij(3)**2
          rmin=min(rmin,sqrt(r2))
          if (r2 < cutoff**2) then
           pair=pair+1 
           f=f+4d0*eps*(sgm12/r2**6-sgm6/r2**3)
          endif
        enddo
      enddo

c-----force
      do i=1,n
        dfdx(1,i)=0d0
        dfdx(2,i)=0d0
        dfdx(3,i)=0d0

        do j=1,n
          if(j.ne.i)then
            rij(1)=x(1,i)-x(1,j)
            rij(2)=x(2,i)-x(2,j)
            rij(3)=x(3,i)-x(3,j)
            rij(1)=rij(1)-h(1,1)*dnint(rij(1)/h(1,1))
            rij(2)=rij(2)-h(2,2)*dnint(rij(2)/h(2,2))
            rij(3)=rij(3)-h(3,3)*dnint(rij(3)/h(3,3))
            r2=rij(1)**2+rij(2)**2+rij(3)**2
            if (r2 < cutoff**2) then
             factor=4d0*eps*
     &         (-12d0*sgm12/r2**7+6d0*sgm6/r2**4)

             dfdx(1,i)=dfdx(1,i)+factor*rij(1)
             dfdx(2,i)=dfdx(2,i)+factor*rij(2)
             dfdx(3,i)=dfdx(3,i)+factor*rij(3)
             do m=1,3
              do l=1,3
              virial(m,l)=virial(m,l)-factor*rij(m)*rij(l)
              enddo
             enddo
            endif
          endif
        enddo
      enddo
      virial=virial*0.5d0

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

! call press      
      subroutine press(v,virial,n,vol,p,ekin)
      implicit none
      real*8 v(3,n),virial(3,3),p(3,3)
      real*8 ekin(3,3),vol
      real*8 amass
      integer i,m,l,n
      parameter(amass=40d0*1836d0)
      ekin=0d0
      do i=1,n
       do m=1,3
       do l=1,3
        ekin(m,l)=ekin(m,l)+amass*v(m,i)*v(l,i)
       enddo
       enddo
      enddo
      p=(ekin+virial)/vol
      end subroutine

!      
