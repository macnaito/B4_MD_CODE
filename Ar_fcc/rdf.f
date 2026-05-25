!boxは手書き

      implicit none
      integer n,i,j,k,nmax,itemp,nbin
      parameter(nmax=30000)
      parameter (nbin=200)
      real*8 x(nmax),y(nmax),z(nmax),tempK(nmax)
      real*8 dx,dy,dz,box,dr,r,rmax,rho,g(nbin)
      character*2 lsp(nmax)


      open(10,file='out001.xyz')
       read(10,*)n
       read(10,*) 
       do i=1,n
        read(10,*)lsp(i),itemp,
     &   x(i),y(i),z(i),tempK(i)
        
       enddo

      box=29.14468d0
      rho=dble(n)/box**3
      rmax=box/2.d0
      dr=rmax/dble(nbin)
      g=0.d0

      do i=1,n-1
        do j=i+1,n
            dx=x(i)-x(j)
            dy=y(i)-y(j)
            dz=z(i)-z(j)
            dx=dx-box*dnint(dx/box)
            dy=dy-box*dnint(dy/box)
            dz=dz-box*dnint(dz/box)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            if(r<rmax) then
                k=int(r/dr)+1
                g(k)=g(k)+2.d0
            endif    
        enddo   
      enddo 

      open(11,file='rdf.dat')
        do i=1,nbin
            r=(i-0.5d0)*dr
            g(i)=g(i)/(4.d0*3.1415926*rho*n*dr*r**2)
            write(11,*)r,g(i)
        enddo
        close(11)
        close(10)
        end





